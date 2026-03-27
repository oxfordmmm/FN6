use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;
use std::{ffi::OsStr, fmt};

use flate2::read::MultiGzDecoder;

use rkyv::{Archive, Deserialize, Serialize};

#[derive(Serialize, Deserialize, Archive)]
#[repr(C)]
pub struct Sample {
    pub header: SampleHeader,
    pub name: String,
    pub is_qc_passed: bool,
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

#[derive(Debug, Serialize, Deserialize, Archive)]
#[repr(C)]
pub struct SampleHeader {
    pub reference_hash: String,
    pub mask_hash: String,
    pub version: String,
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
    for line in reader.lines() {
        if let Ok(ip) = line {
            if !ip.is_empty() {
                if ip.starts_with(">") || ip.starts_with(";") {
                    continue;
                }
                reference.push_str(&ip.trim());
            }
        }
    }

    reference
}

pub fn parse_mask<'a>(filepath: &'a Path, mask: &'a mut Vec<usize>) -> &'a [usize] {
    let reader = get_reader(filepath.to_str().unwrap());
    for line in reader.lines() {
        if let Ok(ip) = line {
            if !ip.is_empty() {
                mask.push(
                    ip.trim()
                        .parse::<usize>()
                        .expect("Failed to parse mask value"),
                );
            }
        }
    }

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
        for line in reader.lines() {
            if let Ok(ip) = line {
                if !ip.is_empty() {
                    if name.is_empty() && ip.starts_with(">") {
                        name = ip
                            .to_string()
                            .split("|")
                            .last()
                            .unwrap()
                            .trim()
                            .replace(">", "");
                        continue;
                    } else if ip.starts_with(">") || ip.starts_with(";") {
                        continue;
                    } else {
                        for ch in ip.chars() {
                            if !(ch as u8 == reference.as_bytes()[char_counter]
                                || mask.binary_search(&char_counter).is_ok())
                            {
                                match ch {
                                    'a' => a.push(char_counter),
                                    'A' => a.push(char_counter),
                                    'c' => c.push(char_counter),
                                    'C' => c.push(char_counter),
                                    'g' => g.push(char_counter),
                                    'G' => g.push(char_counter),
                                    't' => t.push(char_counter),
                                    'T' => t.push(char_counter),
                                    _ => n.push(char_counter),
                                }
                            }
                            char_counter += 1;
                        }
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
            a,
            c,
            g,
            t,
            n,
        }
    }
}

pub fn distance(sample1: &Sample, sample2: &Sample, cutoff: usize) -> Option<usize> {
    if !sample1.is_qc_passed || !sample2.is_qc_passed {
        eprintln!(
            "Neglecting {} and {} because of QC failure",
            sample1.name, sample2.name
        );
        return None;
    }
    let distance = dist(
        sample1.a.clone(),
        sample1.n.clone(),
        sample2.a.clone(),
        sample2.n.clone(),
        cutoff,
    ) + dist(
        sample1.c.clone(),
        sample1.n.clone(),
        sample2.c.clone(),
        sample2.n.clone(),
        cutoff,
    ) + dist(
        sample1.t.clone(),
        sample1.n.clone(),
        sample2.t.clone(),
        sample2.n.clone(),
        cutoff,
    ) + dist(
        sample1.g.clone(),
        sample1.n.clone(),
        sample2.g.clone(),
        sample2.n.clone(),
        cutoff,
    );
    if distance > cutoff {
        return None;
    }

    Some(distance)
}

pub fn dist(
    this_x: Vec<usize>,
    this_n: Vec<usize>,
    sample_x: Vec<usize>,
    sample_n: Vec<usize>,
    cutoff: usize,
) -> usize {
    let mut distance = 0;
    for elem in this_x.iter() {
        if distance > cutoff {
            return distance;
        }
        if sample_x.binary_search(&elem).is_err() {
            // Not in sample
            if sample_n.binary_search(&elem).is_err() {
                // Not an n either
                distance += 1;
            }
        }
    }
    for elem in sample_x.iter() {
        if distance > cutoff {
            return distance;
        }
        if this_x.binary_search(&elem).is_err() {
            // Not in sample
            if this_n.binary_search(&elem).is_err() {
                // Not an n either
                distance += 1;
            }
        }
    }
    distance
}
