//! Module for handling the parsing of the reference, mask and sample FASTA files, as well as the distance calculation between samples.
//!
//! The main struct in this module is the `Sample` struct, which contains the positions of the A, C, G, T and N bases in the sample that differ from the reference and are not masked, as well as some metadata about the sample in the header.
//! Also contains functions for calculating distances between `Sample` structs, or their `rkyv`'ed `ArchivedSample` counterparts, which are used to speed up distance calculations by avoiding the overhead of deserialization.
#![doc = include_str!("../README.md")]
use std::collections::HashSet;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;
use std::{ffi::OsStr, fmt};

use pyo3::prelude::*;

use flate2::read::MultiGzDecoder;

use rkyv::{Archive, Deserialize, Serialize};

/// Header struct to store some metadata about the save
#[pyclass(get_all, set_all, from_py_object, str)]
#[derive(Debug, Serialize, Deserialize, Archive, Clone)]
#[repr(C)]
pub struct SampleHeader {
    /// SHA256 hash of the reference genome
    pub reference_hash: String,

    /// SHA256 hash of the mask
    pub mask_hash: String,

    /// Internal version string to allow for future changes to the save format
    pub version: String,
}

/// Struct to store the compressed sample data.
/// The `a`, `c`, `g`, `t` and `n` fields store the positions of the respective bases in the sample that differ from the reference and are not masked.
/// A Sample is a product of a sample's FASTA file, the reference genome which this is in respect to, and the mask.
#[pyclass(get_all, set_all, str)]
#[derive(Serialize, Deserialize, Archive)]
#[repr(C)]
pub struct Sample {
    /// Header containing metadata about the save
    pub header: SampleHeader,

    /// Name of the sample. Used for reporting what distances are between.
    /// This is either derived from the FASTA header or passed as an argument when creating the sample, and has no effect on the distance calculations.
    pub name: String,

    /// Whether the sample passed QC.
    /// Current QC requires >= 80% ACGT in a sample
    pub is_qc_passed: bool,

    /// Whether this is a legacy save, and so whether the header should be expected to be empty.
    pub is_fn5: bool,

    /// Positions where this sample is A and the reference is not A, and that are not masked.
    pub a: Vec<usize>,

    /// Positions where this sample is C and the reference is not C, and that are not masked.
    pub c: Vec<usize>,

    /// Positions where this sample is T and the reference is not T, and that are not masked.
    pub t: Vec<usize>,

    /// Positions where this sample is G and the reference is not G, and that are not masked.
    pub g: Vec<usize>,

    /// Positions where this sample is non-ACGT.
    pub n: Vec<usize>,
}

impl fmt::Display for Sample {
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

impl fmt::Display for SampleHeader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SampleHeader")
            .field("reference_hash", &self.reference_hash)
            .field("mask_hash", &self.mask_hash)
            .field("version", &self.version)
            .finish()
    }
}

/// Read normal or compressed files seamlessly.
/// Uses the presence of a `gz` extension to choose between the two:
/// <https://users.rust-lang.org/t/write-to-normal-or-gzip-file-transparently/35561>
///
/// # Arguments
/// - `filename` - Path to the file
///
/// # Returns
/// - BufRead object to read the file
fn get_reader(filename: &str) -> Box<dyn BufRead> {
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

/// Parse a given FASTA to extract the string of the genome.
///
/// # Arguments
/// - `filepath` - Path to the FASTA file
///
/// # Returns
/// - String of the genome, with all headers, comments, and newlines removed
///
/// # Panics
/// - If the FASTA file contains more than one header, as this would indicate that it contains more than one contig, which is not supported by this tool
pub fn parse_reference(filepath: &Path) -> String {
    let reader = get_reader(filepath.to_str().unwrap());
    let mut reference = String::new();
    let mut header_count = 0;
    for line in reader.lines().map_while(Result::ok) {
        if !line.is_empty() {
            if line.starts_with(">") {
                header_count += 1;
                continue;
            }
            if line.starts_with(";") {
                continue;
            }
            reference.push_str(line.trim());
        }
    }
    if header_count != 1 {
        panic!(
            "Reference file {} contains more than one header. Only FASTA files of a single contig are allowed.",
            filepath.display()
        );
    }

    reference
}

/// Parse a given mask file to extract the positions to be masked.
/// Positions should be given as a line separated file, with 0-based positions.
///
/// # Arguments
/// - `filepath` - Path to the mask file
///
/// # Returns
/// - Vector of positions to be masked, sorted in ascending order
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
    /// Instanciate a new sample from a FASTA file, given the reference and mask to compare against.
    /// The sample will be parsed and the positions of the A, C, G, T and N bases that differ from the reference and are not masked will be stored in the respective fields of the struct.
    ///
    /// # Arguments
    /// - `filepath` - Path to the sample FASTA file
    /// - `reference` - String of the reference genome to compare against
    /// - `mask` - Vector of positions to be masked (i.e., ignored) during the analysis
    /// - `mask_hash` - SHA256 hash of the mask, for storing in the header
    /// - `reference_hash` - SHA256 hash of the reference, for storing in the header
    ///
    /// # Returns
    /// - Sample struct containing the parsed data from the FASTA file and the metadata in the header
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
        // Technically these should be sorted for free given everything is incrementing
        // But double check for sanity as this should only be done once
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

/// Parse a given FN5 file to extract the sample data and convert it to the new format.
/// Returns the bytes format of the ArchivedSample struct, which can be directly compared to the bytes format of samples created with the new method.
/// This extracts the sample name from the filepath, matching the FN5 format.
///
/// # Arguments
/// - `filepath` - Path to the FN5 file
///
/// # Returns
/// - Vector of bytes representing the ArchivedSample struct containing the parsed data from the FN5 file
///
/// # Panics
/// - If the FN5 file is not in the correct format:
///     - If its length is not a multiple of 4, as it should consist of 32-bit integers
///     - If the number of positions to read for any of the A, C, G, T or N fields does not match the actual number of positions provided for that field
pub fn from_fn5(filepath: &Path) -> Vec<u8> {
    let name = filepath.file_stem().unwrap().to_str().unwrap().to_string();
    let mut results = [Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()];
    // Read full file into RAM before trying to do anything
    let bytes = std::fs::read(filepath).unwrap();
    if !bytes.len().is_multiple_of(4) {
        panic!("File {} is not in the correct format", filepath.display());
    }

    /* FN5 format is as follows:

    - 4 bytes: number of positions to read for A field
    - 4 bytes for each position in A field
    - 4 bytes: number of positions to read for C field
    - 4 bytes for each position in C field
    - 4 bytes: number of positions to read for T field
    - 4 bytes for each position in T field
    - 4 bytes: number of positions to read for G field
    - 4 bytes for each position in G field
    - 4 bytes: number of positions to read for N field
    - 4 bytes for each position in N field

    Each 4 byte chunk is a little-endian 32-bit integer
     */

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

    if chunks_to_read != 0 {
        panic!("File {} is not in the correct format", filepath.display());
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

/// Calculate the distance between two samples, defined as the number of positions where they differ and that are not masked.
/// This is done by comparing the positions of the A, C, G, T and N fields of the two samples, and counting the number of positions that are different between them.
/// The distance is returned as an Option, where None indicates that the distance is above the given cutoff, and Some(distance) indicates that the distance is below the cutoff.
/// The cutoff is used to speed up the distance calculation by allowing for early termination if the distance is already above the cutoff, as this is often the case when comparing samples that are not closely related.
///
/// # Arguments
/// - `sample1` - First sample to compare
/// - `sample2` - Second sample to compare
/// - `cutoff` - Maximum distance to calculate before returning None
///
/// # Returns
/// - `Option<usize>` representing the distance between the two samples, or None if the distance is above the cutoff
///
/// # Panics
/// - If the two samples are in respect to different masks or references, as this would make the distance calculation meaningless. Note that this is skipped if one of the saves is from FN5
/// - If the two samples are in respect to different versions, as this would indicate that they are not comparable. Note that this is skipped if one of the saves is from FN5
#[pyfunction(signature = (sample1, sample2, cutoff = 20))]
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

/// Calculate the distance between two samples, defined as the number of positions where they differ and that are not masked.
/// Notably this differs from `distance` in that it takes ArchivedSamples as input, which are generated by `rkyv` as a way to avoid the overhead of deserialization.
/// This is done by comparing the positions of the A, C, G, T and N fields of the two samples, and counting the number of positions that are different between them.
/// The distance is returned as an Option, where None indicates that the distance is above the given cutoff, and Some(distance) indicates that the distance is below the cutoff.
/// The cutoff is used to speed up the distance calculation by allowing for early termination if the distance is already above the cutoff, as this is often the case when comparing samples that are not closely related.
///
/// # Arguments
/// - `sample1` - First sample to compare
/// - `sample2` - Second sample to compare
/// - `cutoff` - Maximum distance to calculate before returning None
///
/// # Returns
/// - `Option<usize>` representing the distance between the two samples, or None if the distance is above the cutoff
///
/// # Panics
/// - If the two samples are in respect to different masks or references, as this would make the distance calculation meaningless. Note that this is skipped if one of the saves is from FN5
/// - If the two samples are in respect to different versions, as this would indicate that they are not comparable. Note that this is skipped if one of the saves is from FN5
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

/// Helper function to calculate the distance between two samples for a given base, defined as the number of positions where they differ and that are not masked.
///
/// # Arguments
/// - `this_x` - Positions of the base in the first sample
/// - `this_n` - Positions of the N bases in the first sample
/// - `sample_x` - Positions of the base in the second sample
/// - `sample_n` - Positions of the N bases in the second sample
/// - `cutoff` - Maximum distance to calculate before returning
/// - `distances` - HashSet to store the positions where the two samples differ for the given base, which is used to calculate the distance between the two samples. Note that this is required otherwise SNPs can be double counted if the same position is different between the two samples for multiple bases (e.g., A in one sample and C in the other, which would be counted as a difference for both the A and C fields if the positions were not stored in a HashSet to prevent double counting).
fn dist<T: std::cmp::Ord + std::hash::Hash + Copy>(
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;

    use super::*;

    #[test]
    fn test_dist() {
        let mut distances = HashSet::new();
        // Basic dists
        dist(&[1, 2, 3], &[4], &[2, 3, 4], &[], 10, &mut distances);
        assert_eq!(distances.len(), 1);
        assert!(distances.contains(&1));

        distances.clear();

        dist(
            &[1, 2, 3, 4, 5, 6],
            &[4],
            &[2, 3, 4],
            &[5],
            10,
            &mut distances,
        );
        assert_eq!(distances.len(), 2);
        assert!(distances.contains(&1));
        assert!(distances.contains(&6));

        distances.clear();

        // Check that cutoff works
        // There should be 2 SNPs here, but with a cutoff of 0, it should return at 1 SNP (to distinguish between dist at cutoff and above cutoff)
        dist(
            &[1, 2, 3, 4, 5, 6],
            &[4],
            &[2, 3, 4],
            &[5],
            0,
            &mut distances,
        );
        assert_eq!(distances.len(), 1);

        distances.clear();

        // Similarly but with 3 SNPs and cutoff 1
        dist(
            &[1, 2, 3, 4, 5, 6, 7],
            &[4],
            &[2, 3, 4],
            &[5],
            1,
            &mut distances,
        );
        assert_eq!(distances.len(), 2);
    }
}
