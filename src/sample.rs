use std::collections::HashSet;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;
use std::{ffi::OsStr, fmt};

use flate2::read::MultiGzDecoder;

use rkyv::{Archive, Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, Archive)]
#[repr(C)]
pub struct SampleHeader {
    pub reference_hash: String,
    pub mask_hash: String,
    pub version: String,
}
#[derive(Serialize, Deserialize, Archive)]
#[repr(C)]
pub struct Sample {
    pub header: SampleHeader,
    pub name: String,
    pub is_qc_passed: bool,
    pub is_fn5: bool,
    pub a: Vec<usize>,
    pub c: Vec<usize>,
    pub t: Vec<usize>,
    pub g: Vec<usize>,
    pub n: Vec<usize>,
}

impl fmt::Debug for Sample {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Sample")
            .field("header", &self.header)
            .field("name", &self.name)
            .field("is_qc_passed", &self.is_qc_passed)
            .field("a", &self.a.len())
            .field("c", &self.c.len())
            .field("g", &self.t.len())
            .field("t", &self.g.len())
            .field("n", &self.n.len())
            .finish()
    }
}

/// Read normal or compressed files seamlessly
/// Uses the presence of a `gz` extension to choose between the two
/// https://users.rust-lang.org/t/write-to-normal-or-gzip-file-transparently/35561
///
/// # Arguments
/// - `filename` - Path to the file
///
/// # Returns
/// - BufRead object to read the file
pub fn get_reader(filename: &str) -> Box<dyn BufRead> {
    let path = Path::new(filename);
    let file = match File::open(path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why),
        Ok(file) => file,
    };

    if path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufReader::with_capacity(
            128 * 1024,
            MultiGzDecoder::new(file),
        ))
    } else {
        Box::new(BufReader::with_capacity(128 * 1024, file))
    }
}

pub fn parse_reference(filepath: &Path) -> String {
    let reader = get_reader(filepath.to_str().unwrap());
    let mut reference = String::new();
    for line in reader.lines().map_while(Result::ok) {
        if !line.is_empty() {
            if line.starts_with(">") || line.starts_with(";") {
                continue;
            }
            reference.push_str(line.trim());
        }
    }

    reference
}

pub fn parse_mask(filepath: &Path) -> Vec<usize> {
    let reader = get_reader(filepath.to_str().unwrap());
    let mut mask = Vec::new();
    for line in reader.lines().map_while(Result::ok) {
        if !line.is_empty() {
            mask.push(
                line.trim()
                    .parse::<usize>()
                    .expect("Failed to parse mask value"),
            );
        }
    }
    mask.sort();

    mask
}

impl Sample {
    pub fn new(
        filepath: &Path,
        reference: &str,
        mask: &[usize],
        mask_hash: &str,
        reference_hash: &str,
    ) -> Sample {
        let reader = get_reader(filepath.to_str().unwrap());
        let mut name = String::new();
        let mut a = Vec::new();
        let mut c = Vec::new();
        let mut g = Vec::new();
        let mut t = Vec::new();
        let mut n = Vec::new();

        let mut char_counter = 0;
        for line in reader.lines().map_while(Result::ok) {
            if !line.is_empty() {
                if name.is_empty() && line.starts_with(">") {
                    // This is a bit hacky but allows for some variation in the header format
                    // `>chr1|chr1.fa`, `>chr1 description` and `>chr1` will all work, and the sample name will be `chr1` in both cases
                    name = line
                        .split_whitespace()
                        .collect::<Vec<&str>>()
                        .first()
                        .unwrap()
                        .split("|")
                        .last()
                        .unwrap()
                        .trim()
                        .replace(">", "");
                    continue;
                } else if line.starts_with(">") || line.starts_with(";") {
                    continue;
                } else {
                    for ch in line.chars() {
                        if !(ch as u8 == reference.as_bytes()[char_counter]
                            || mask.binary_search(&char_counter).is_ok())
                        {
                            match ch {
                                'A' | 'a' => a.push(char_counter),
                                'C' | 'c' => c.push(char_counter),
                                'G' | 'g' => g.push(char_counter),
                                'T' | 't' => t.push(char_counter),
                                _ => n.push(char_counter),
                            }
                        }
                        char_counter += 1;
                    }
                }
            }
        }
        a.sort_unstable();
        c.sort_unstable();
        g.sort_unstable();
        t.sort_unstable();
        n.sort_unstable();

        Sample {
            // Metadata about the sample's save
            // Allows for checking compatibility of samples for comparison
            header: SampleHeader {
                reference_hash: reference_hash.to_string(),
                mask_hash: mask_hash.to_string(),
                version: "1.0".to_string(),
            },
            name,
            is_qc_passed: (n.len() as f64 / reference.len() as f64) < 0.2,
            is_fn5: false,
            a,
            c,
            g,
            t,
            n,
        }
    }
}

pub fn from_fn5(filepath: &Path) -> Vec<u8> {
    let name = filepath.file_stem().unwrap().to_str().unwrap().to_string();
    let mut results = [Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()];
    // Read full file into RAM before trying to do anything
    let bytes = std::fs::read(filepath).unwrap();
    if !bytes.len().is_multiple_of(4) {
        panic!("File {} is not in the correct format", filepath.display());
    }

    let mut idx: i32 = -1;
    let mut chunks_to_read = 0;
    for chunk in bytes.chunks(4) {
        let val = i32::from_le_bytes(chunk.try_into().unwrap());
        if chunks_to_read == 0 {
            chunks_to_read = val;
            idx += 1;
        } else {
            results[idx as usize].push(val as usize);
            chunks_to_read -= 1;
        }
    }

    let sample = Sample {
        header: SampleHeader {
            reference_hash: String::new(),
            mask_hash: String::new(),
            version: String::new(),
        },
        name,
        is_qc_passed: true,
        is_fn5: true,
        a: results[0].clone(),
        c: results[1].clone(),
        g: results[2].clone(),
        t: results[3].clone(),
        n: results[4].clone(),
    };
    rkyv::to_bytes::<rkyv::rancor::Error>(&sample)
        .unwrap()
        .to_vec()
}

pub fn distance(sample1: &Sample, sample2: &Sample, cutoff: usize) -> Option<usize> {
    if !sample1.is_qc_passed || !sample2.is_qc_passed {
        eprintln!(
            "Neglecting {} and {} because of QC failure",
            sample1.name, sample2.name
        );
        return None;
    }
    // Be flexible with whether a sample requires a header
    // This is to enable comparison of samples from FN5 which don't have a header
    if !sample1.is_fn5 && !sample2.is_fn5 {
        let header1 = &sample1.header;
        let header2 = &sample2.header;
        if header1.mask_hash != header2.mask_hash {
            panic!(
                "Cannot compare {} and {} because of different masks",
                sample1.name, sample2.name
            );
        }
        if header1.reference_hash != header2.reference_hash {
            panic!(
                "Cannot compare {} and {} because of different references",
                sample1.name, sample2.name
            );
        }
        if header1.version != header2.version {
            panic!(
                "Cannot compare {} and {} because of different versions",
                sample1.name, sample2.name
            );
        }
    }
    let mut distances = HashSet::new();
    dist(
        &sample1.a,
        &sample1.n,
        &sample2.a,
        &sample2.n,
        cutoff,
        &mut distances,
    );
    dist(
        &sample1.c,
        &sample1.n,
        &sample2.c,
        &sample2.n,
        cutoff,
        &mut distances,
    );
    dist(
        &sample1.t,
        &sample1.n,
        &sample2.t,
        &sample2.n,
        cutoff,
        &mut distances,
    );
    dist(
        &sample1.g,
        &sample1.n,
        &sample2.g,
        &sample2.n,
        cutoff,
        &mut distances,
    );
    if distances.len() > cutoff {
        return None;
    }

    Some(distances.len())
}

pub fn arch_distance(
    sample1: &ArchivedSample,
    sample2: &ArchivedSample,
    cutoff: usize,
) -> Option<usize> {
    if !sample1.is_qc_passed || !sample2.is_qc_passed {
        eprintln!(
            "Neglecting {} and {} because of QC failure",
            sample1.name, sample2.name
        );
        return None;
    }
    // Be flexible with whether a sample requires a header
    // This is to enable comparison of samples from FN5 which don't have a header
    if !sample1.is_fn5 && !sample2.is_fn5 {
        let header1 = &sample1.header;
        let header2 = &sample2.header;
        if header1.mask_hash != header2.mask_hash {
            panic!(
                "Cannot compare {} and {} because of different masks",
                sample1.name, sample2.name
            );
        }
        if header1.reference_hash != header2.reference_hash {
            panic!(
                "Cannot compare {} and {} because of different references",
                sample1.name, sample2.name
            );
        }
        if header1.version != header2.version {
            panic!(
                "Cannot compare {} and {} because of different versions",
                sample1.name, sample2.name
            );
        }
    }
    let mut distances = HashSet::new();
    dist(
        &sample1.a,
        &sample1.n,
        &sample2.a,
        &sample2.n,
        cutoff,
        &mut distances,
    );
    dist(
        &sample1.c,
        &sample1.n,
        &sample2.c,
        &sample2.n,
        cutoff,
        &mut distances,
    );
    dist(
        &sample1.t,
        &sample1.n,
        &sample2.t,
        &sample2.n,
        cutoff,
        &mut distances,
    );
    dist(
        &sample1.g,
        &sample1.n,
        &sample2.g,
        &sample2.n,
        cutoff,
        &mut distances,
    );
    if distances.len() > cutoff {
        return None;
    }

    Some(distances.len())
}

pub fn dist<T: std::cmp::Ord + std::hash::Hash + Copy>(
    this_x: &[T],
    this_n: &[T],
    sample_x: &[T],
    sample_n: &[T],
    cutoff: usize,
    distances: &mut HashSet<T>,
) {
    for elem in this_x.iter() {
        if distances.len() > cutoff {
            return;
        }
        if sample_x.binary_search(elem).is_err() {
            // Not in sample
            if sample_n.binary_search(elem).is_err() {
                // Not an n either
                distances.insert(*elem);
            }
        }
    }
    for elem in sample_x.iter() {
        if distances.len() > cutoff {
            return;
        }
        if this_x.binary_search(elem).is_err() {
            // Not in sample
            if this_n.binary_search(elem).is_err() {
                // Not an n either
                distances.insert(*elem);
            }
        }
    }
}
