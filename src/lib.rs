use std::{
    path::{Path, PathBuf},
    sync::Mutex,
};

use crate::sample::ArchivedSample;
use rayon::prelude::*;

pub mod sample;

/// Seemlessly load either a new fasta or an existing save
/// This is used by the CLI to allow users to either load existing .fn6 files or create new ones on the fly. If a .fn6 file is provided, it will be loaded and deserialized. If a .fn5 file is provided, it will be loaded and deserialized. If a FASTA file is provided, it will be processed and compressed into a Sample struct.
pub fn load_save(
    filepath: &Path,
    reference: &str,
    mask: &[usize],
    mask_hash: &str,
    reference_hash: &str,
) -> sample::Sample {
    if filepath.to_str().unwrap().to_owned().ends_with(".fn6") {
        let bytes = std::fs::read(filepath).unwrap();
        let arch_sample = rkyv::access::<ArchivedSample, rkyv::rancor::Error>(&bytes[..]).unwrap();
        // Deserialising is relatively expensive - using this leads to double RAM usage + slower runtime.
        return rkyv::deserialize::<sample::Sample, rkyv::rancor::Error>(arch_sample).unwrap();
    }
    if filepath.to_str().unwrap().to_owned().ends_with(".fn5") {
        let bytes = sample::from_fn5(filepath);
        let arch_sample = rkyv::access::<ArchivedSample, rkyv::rancor::Error>(&bytes[..]).unwrap();
        // Deserialising is relatively expensive - using this leads to double RAM usage + slower runtime.
        return rkyv::deserialize::<sample::Sample, rkyv::rancor::Error>(arch_sample).unwrap();
    }
    sample::Sample::new(filepath, reference, mask, mask_hash, reference_hash)
}

/// Reference compress a sample (parse into `sample::Sample`) and save it to disk if it passes QC.
///
/// # Arguments
/// - `filepath`: Path to the sample genome FASTA file
/// - `reference`: The reference genome sequence as a string. This is used for compression and distance calculations.
/// - `mask`: A list of positions in the reference genome that should be masked (i.e., ignored) during the analysis. The positions are 0-based.
/// - `mask_hash`: A hash of the mask file. This is used for QC to ensure that the same mask is used for all samples.
/// - `reference_hash`: A hash of the reference genome. This is used for QC to ensure that the same reference is used for all samples.
///
/// # Optional Arguments
/// - `id`: An optional ID for this sample. If not provided, the ID will be derived from the header of the FASTA file.
/// - `output_path`: An optional output path for the .fn6 file. If not provided, the .fn6 file will be saved in the same directory as the sample FASTA file with the same name but with a .fn6 extension.
///
/// # Returns
/// The `sample::Sample` struct representing the compressed sample. This struct contains the compressed representation of the sample genome, as well as metadata about the reference and mask used for compression, and whether the sample passed QC. If used in downstream analysis, it is the user's responsibility to check that the sample passed QC by checking the `is_qc_passed` field of the returned `sample::Sample` struct.
pub fn reference_compress(
    filepath: &Path,
    reference: &str,
    mask: &[usize],
    mask_hash: &str,
    reference_hash: &str,
    id: Option<String>,
    output_path: Option<PathBuf>,
) -> sample::Sample {
    let mut s = sample::Sample::new(filepath, reference, mask, mask_hash, reference_hash);
    if let Some(name) = id {
        s.name = name;
    }
    let output = match output_path {
        Some(path) => path,
        None => {
            let mut path = filepath.to_path_buf();
            path.set_extension("fn6");
            path
        }
    };
    // Only save if qc is passed
    if s.is_qc_passed {
        let serialized = rkyv::to_bytes::<rkyv::rancor::Error>(&s).unwrap();
        std::fs::write(output, serialized).unwrap();
    }

    s
}

/// Load a set of .fn6 or .fn5 files into memory as byte vectors. This uses multithreading for performance.
/// The byte vectors returned correspond to the `rkyv` serialized `sample::Sample` structs.
///
/// # Arguments
/// - `filepaths`: A vector of paths to .fn6 or .fn5 files to load. The file type is determined by the file extension. .fn6 files are expected to be `rkyv` serialized `sample::Sample` structs, while .fn5 files are expected to be in the old format and will be converted to `sample::Sample` structs and then serialized using `rkyv`.
///
/// # Returns
/// A vector of byte vectors, where each byte vector is the `rkyv` serialized representation of a `sample::Sample` struct. The order of the byte vectors corresponds to the order of the input filepaths.
fn load_arch_saves(filepaths: Vec<PathBuf>) -> Vec<Vec<u8>> {
    filepaths
        .par_iter()
        .map(|sample_path| {
            if sample_path.extension().and_then(|s| s.to_str()) == Some("fn6") {
                std::fs::read(sample_path).unwrap()
            } else if sample_path.extension().and_then(|s| s.to_str()) == Some("fn5") {
                sample::from_fn5(sample_path)
            } else {
                panic!("Unsupported file type: {:?}", sample_path);
            }
        })
        .collect()
}

/// Given a set of comparisons to do, compute the distances and print them to stdout as they are computed. This uses multithreading for performance, and a mutex to ensure that the output is not interleaved.
///
/// # Arguments
/// - `comparisons`: A vector of tuples, where each tuple contains two byte vectors. Each byte vector is the `rkyv` serialized representation of a `sample::Sample` struct. The order of the byte vectors corresponds to the order of the input filepaths used to load the samples.
/// - `cutoff`: The SNP threshold for distance calculation. If the distance between two samples exceeds this threshold, the distance will not be reported. This is used to speed up distance calculations by allowing for early termination when the distance is large.
///
/// # Output
/// The function prints the distances to stdout as they are computed. Each line of output corresponds to a pairwise comparison and is formatted as "sample1_name sample2_name distance". If the distance exceeds the cutoff, no distance is written.
fn get_distances(comparisons: Vec<(&Vec<u8>, &Vec<u8>)>, cutoff: usize) {
    let distances: Mutex<Vec<(String, String, usize)>> = Mutex::new(Vec::new());

    let _ = comparisons
        .par_iter()
        .map(|(sample1, sample2)| {
            let sample1 =
                rkyv::access::<ArchivedSample, rkyv::rancor::Error>(&sample1[..]).unwrap();
            let sample2 =
                rkyv::access::<ArchivedSample, rkyv::rancor::Error>(&sample2[..]).unwrap();
            if sample1.name == sample2.name {
                return;
            }
            let dist = sample::arch_distance(sample1, sample2, cutoff);
            if let Some(d) = dist {
                let mut dist_lock = distances.lock().unwrap();
                dist_lock.push((sample1.name.to_string(), sample2.name.to_string(), d));
                if dist_lock.len() == 1000 {
                    for (name1, name2, d) in dist_lock.iter() {
                        println!("{} {} {}", name1, name2, d);
                    }
                    dist_lock.clear();
                }
            }
        })
        .collect::<Vec<()>>();

    // Catch the last bit of distances
    let dist_lock = distances.lock().unwrap();
    for (name1, name2, d) in dist_lock.iter() {
        println!("{} {} {}", name1, name2, d);
    }
}

/// Compute all distances from a vec of genome save file paths. This is the main function for the "compute" command in the CLI. It loads the samples, figures out what comparisons to do, and then calls `get_distances` to compute and print the distances.
///
/// # Arguments
/// - `filepaths`: A vector of paths to .fn6 or .fn5 files to load. The file type is determined by the file extension. .fn6 files are expected to be `rkyv` serialized `sample::Sample` structs, while .fn5 files are expected to be in the old format and will be converted to `sample::Sample` structs and then serialized using `rkyv`.
/// - `cutoff`: The SNP threshold for distance calculation. If the distance between two samples exceeds this threshold, it will not be reported. This is used to speed up distance calculations by allowing for early termination when the distance is large.
pub fn compute(filepaths: Vec<PathBuf>, cutoff: usize) {
    // Load the saves
    let samples = load_arch_saves(filepaths);

    // Figure out what comparisons we need to do
    let mut comparisons: Vec<(&Vec<u8>, &Vec<u8>)> = Vec::new();
    for (idx, sample1) in samples.iter().enumerate() {
        for sample2 in samples.iter().skip(idx + 1) {
            comparisons.push((sample1, sample2));
        }
    }

    get_distances(comparisons, cutoff);
}

/// Compute the distances required to add new samples to an existing set. This is the main function for the "add-samples" command in the CLI. It loads the existing and new samples, figures out what comparisons to do, and then calls `get_distances` to compute and print the distances.
///
/// # Arguments
/// - `existing`: A vector of paths to existing .fn6 or .fn5 files to load. The file type is determined by the file extension. .fn6 files are expected to be `rkyv` serialized `sample::Sample` structs, while .fn5 files are expected to be in the old format and will be converted to `sample::Sample` structs and then serialized using `rkyv`.
/// - `new_samples`: A vector of paths to new .fn6 or .fn5 files to load. The file type is determined by the file extension. .fn6 files are expected to be `rkyv` serialized `sample::Sample` structs, while .fn5 files are expected to be in the old format and will be converted to `sample::Sample` structs and then serialized using `rkyv`.
/// - `cutoff`: The SNP threshold for distance calculation. If the distance between two samples exceeds this threshold, it will not be reported. This is used to speed up distance calculations by allowing for early termination when the distance is large.
pub fn add_samples(existing: Vec<PathBuf>, new_samples: Vec<PathBuf>, cutoff: usize) {
    let existing_samples = load_arch_saves(existing);
    let new_samples = load_arch_saves(new_samples);

    let mut comparisons: Vec<(&Vec<u8>, &Vec<u8>)> = Vec::new();
    // Compare each existing sample to each new sample
    for sample1 in existing_samples.iter() {
        for sample2 in new_samples.iter() {
            comparisons.push((sample1, sample2));
        }
    }
    // And each new sample to each other
    for (idx, sample1) in new_samples.iter().enumerate() {
        for sample2 in new_samples.iter().skip(idx + 1) {
            comparisons.push((sample1, sample2));
        }
    }

    get_distances(comparisons, cutoff);
}
