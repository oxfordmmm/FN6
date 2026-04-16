//! Fast, efficient SNP distance calculation from disk.
//!
//! Approximately 10x faster than FN5, easier to use and maintain, and adds checking for matching reference and mask in FN6 saves, all while retaining interoperability with FN5 saves.
use std::{
    path::{Path, PathBuf},
    sync::Mutex,
};

use crate::sample::ArchivedSample;
use rayon::prelude::*;

pub mod py_lib;

pub mod sample;

/// Seemlessly load either a new fasta or an existing save.
/// This is used by the CLI to allow users to either load existing .fn6 files or create new ones on the fly. If a .fn6 file is provided, it will be loaded and deserialized. If a .fn5 file is provided, it will be loaded and deserialized. If a FASTA file is provided, it will be processed and compressed into a Sample struct.
/// All saves are then serialised to bytes with rkyv.
pub fn load_save(
    filepath: &Path,
    reference: &str,
    mask: &[usize],
    mask_hash: &str,
    reference_hash: &str,
) -> Vec<u8> {
    if filepath.to_str().unwrap().to_owned().ends_with(".fn6") {
        return std::fs::read(filepath).unwrap();
    }
    if filepath.to_str().unwrap().to_owned().ends_with(".fn5") {
        return sample::from_fn5(filepath);
    }
    let s = sample::Sample::new(filepath, reference, mask, mask_hash, reference_hash);
    rkyv::to_bytes::<rkyv::rancor::Error>(&s).unwrap().to_vec()
}

/// Load a set of files into memory as byte vectors. This uses multithreading for performance.
/// The byte vectors returned correspond to the `rkyv` serialized `sample::Sample` structs.
///
/// # Arguments
/// - `filepaths`: A vector of paths to load. The file type is determined by the file extension. .fn6 files are expected to be `rkyv` serialized `sample::Sample` structs, while .fn5 files are expected to be in the old format and will be converted to `sample::Sample` structs and then serialized using `rkyv`. FASTA files will be processed and compressed into `sample::Sample` structs and then serialized using `rkyv`.
/// - `reference`: The reference genome sequence as a string. This is only required if at least 1 FASTA file is input
/// - `mask`: A list of positions in the reference genome that should be masked (i.e., ignored) during the analysis. The positions are 0-based. This is only required if at least 1 FASTA file is input.
/// - `mask_hash`: A hash of the mask file. This is used for QC to ensure that the same mask is used for all samples. This is only required if at least 1 FASTA file is input.
/// - `reference_hash`: A hash of the reference genome. This is used for QC to ensure that the same reference is used for all samples. This is only required if at least 1 FASTA file is input.
///
/// # Returns
/// A vector of byte vectors, where each byte vector is the `rkyv` serialized representation of a `sample::Sample` struct. The order of the byte vectors corresponds to the order of the input filepaths.
pub fn load_arch_saves(
    filepaths: Vec<PathBuf>,
    reference: &str,
    mask: &[usize],
    mask_hash: &str,
    reference_hash: &str,
) -> Vec<Vec<u8>> {
    filepaths
        .par_iter()
        .map(|sample_path| load_save(sample_path, reference, mask, mask_hash, reference_hash))
        .collect()
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

/// Given a set of comparisons to do, compute the distances and print them to stdout as they are computed. This uses multithreading for performance, and a mutex to ensure that the output is not interleaved.
///
/// # Arguments
/// - `comparisons`: A vector of tuples, where each tuple contains two byte vectors. Each byte vector is the `rkyv` serialized representation of a `sample::Sample` struct. The order of the byte vectors corresponds to the order of the input filepaths used to load the samples.
/// - `cutoff`: The SNP threshold for distance calculation. If the distance between two samples exceeds this threshold, the distance will not be reported. This is used to speed up distance calculations by allowing for early termination when the distance is large.
///
/// # Output
/// The function prints the distances to stdout as they are computed. Each line of output corresponds to a pairwise comparison and is formatted as "sample1_name sample2_name distance". If the distance exceeds the cutoff, no distance is written.
pub fn get_distances(comparisons: Vec<(&Vec<u8>, &Vec<u8>)>, cutoff: usize) {
    let distances: Mutex<Vec<(String, String, usize)>> = Mutex::new(Vec::new());
    distances.lock().unwrap().reserve(1000);

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
/// - `filepaths`: A vector of files to load. The file type is determined by the file extension. .fn6 files are expected to be `rkyv` serialized `sample::Sample` structs, while .fn5 files are expected to be in the old format and will be converted to `sample::Sample` structs and then serialized using `rkyv`. FASTA files will be reference compressed on the fly and then serialized using `rkyv`.
/// - `reference`: The reference genome sequence as a string. This is only required if at least 1 FASTA file is input
/// - `mask`: A list of positions in the reference genome that should be masked (i.e., ignored) during the analysis. The positions are 0-based. This is only required if at least 1 FASTA file is input.
/// - `mask_hash`: A hash of the mask file. This is used for QC to ensure that the same mask is used for all samples. This is only required if at least 1 FASTA file is input.
/// - `reference_hash`: A hash of the reference genome. This is used for QC to ensure that the same reference is used for all samples. This is only required if at least 1 FASTA file is input.
/// - `cutoff`: The SNP threshold for distance calculation. If the distance between two samples exceeds this threshold, it will not be reported. This is used to speed up distance calculations by allowing for early termination when the distance is large.
pub fn compute(
    filepaths: Vec<PathBuf>,
    reference: &str,
    mask: &[usize],
    mask_hash: &str,
    reference_hash: &str,
    cutoff: usize,
    debug: bool,
) {
    // Load the saves
    let start_time = std::time::Instant::now();
    let samples = load_arch_saves(
        filepaths.clone(),
        reference,
        mask,
        mask_hash,
        reference_hash,
    );
    if debug {
        eprintln!(
            "Loaded {} samples in {:.2?}",
            samples.len(),
            start_time.elapsed()
        );
    }

    // Figure out what comparisons we need to do
    let mut n_comps: u64 = 0;
    let mut comparisons: Vec<(&Vec<u8>, &Vec<u8>)> = Vec::new();
    for (idx, sample1) in samples.iter().enumerate() {
        for sample2 in samples.iter().skip(idx + 1) {
            comparisons.push((sample1, sample2));
            if comparisons.len() > 1_000_000 {
                // We're doing a lot of comparisons, so batch them up to avoid excessive RAM usage
                get_distances(comparisons, cutoff);
                comparisons = Vec::new();
            }
            n_comps += 1;
        }
    }

    // Get last distances
    get_distances(comparisons, cutoff);

    if debug {
        eprintln!(
            "Computed {} distances in {:.2?} ({:.2?}) per comparison",
            n_comps,
            start_time.elapsed(),
            start_time.elapsed() / n_comps as u32
        );
    }
}

/// Compute the distances required to add new samples to an existing set. This is the main function for the "add-samples" command in the CLI. It loads the existing and new samples, figures out what comparisons to do, and then calls `get_distances` to compute and print the distances.
///
/// # Arguments
/// - `existing`: A vector of files to load. The file type is determined by the file extension. .fn6 files are expected to be `rkyv` serialized `sample::Sample` structs, while .fn5 files are expected to be in the old format and will be converted to `sample::Sample` structs and then serialized using `rkyv`. FASTA files will be reference compressed on the fly and then serialized using `rkyv`. These are the existing samples that we want to add to.
/// - `new_samples`: A vector of files to load. The file type is determined by the file extension. .fn6 files are expected to be `rkyv` serialized `sample::Sample` structs, while .fn5 files are expected to be in the old format and will be converted to `sample::Sample` structs and then serialized using `rkyv`. FASTA files will be reference compressed on the fly and then serialized using `rkyv`. These are the new samples that we want to add to the existing set.
/// - `reference`: The reference genome sequence as a string. This is only required if at least 1 FASTA file is input
/// - `mask`: A list of positions in the reference genome that should be masked (i.e., ignored) during the analysis. The positions are 0-based. This is only required if at least 1 FASTA file is input.
/// - `mask_hash`: A hash of the mask file. This is used for QC to ensure that the same mask is used for all samples. This is only required if at least 1 FASTA file is input.
/// - `reference_hash`: A hash of the reference genome. This is used for QC to ensure that the same reference is used for all samples. This is only required if at least 1 FASTA file is input.
/// - `cutoff`: The SNP threshold for distance calculation. If the distance between two samples exceeds this threshold, it will not be reported. This is used to speed up distance calculations by allowing for early termination when the distance is large.
#[allow(clippy::too_many_arguments)] // It's not that many and they're all important...
pub fn add_samples(
    existing: Vec<PathBuf>,
    new_samples: Vec<PathBuf>,
    reference: &str,
    mask: &[usize],
    mask_hash: &str,
    reference_hash: &str,
    cutoff: usize,
    debug: bool,
) {
    let start_time = std::time::Instant::now();
    let existing_samples = load_arch_saves(existing, reference, mask, mask_hash, reference_hash);
    let new_samples = load_arch_saves(new_samples, reference, mask, mask_hash, reference_hash);

    let mut comparisons: Vec<(&Vec<u8>, &Vec<u8>)> = Vec::new();
    let mut n_comps: u64 = 0;
    // Compare each existing sample to each new sample
    for sample1 in existing_samples.iter() {
        for sample2 in new_samples.iter() {
            comparisons.push((sample1, sample2));
            if comparisons.len() > 1_000_000 {
                // We're doing a lot of comparisons, so batch them up to avoid excessive RAM usage
                get_distances(comparisons, cutoff);
                comparisons = Vec::new();
            }
            n_comps += 1;
        }
    }
    // And each new sample to each other
    for (idx, sample1) in new_samples.iter().enumerate() {
        for sample2 in new_samples.iter().skip(idx + 1) {
            comparisons.push((sample1, sample2));
            if comparisons.len() > 1_000_000 {
                // We're doing a lot of comparisons, so batch them up to avoid excessive RAM usage
                get_distances(comparisons, cutoff);
                comparisons = Vec::new();
            }
            n_comps += 1;
        }
    }

    get_distances(comparisons, cutoff);

    if debug {
        eprintln!(
            "Computed {} new distances in {:.2?} ({:.2?}) per comparison",
            n_comps,
            start_time.elapsed(),
            start_time.elapsed() / n_comps as u32
        );
    }
}
