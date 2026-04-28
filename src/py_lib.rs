#![cfg(not(tarpaulin_include))]
//! Fast, efficient SNP distance calculation from disk.
//!
//! Approximately 10x faster than FN5, easier to use and maintain, and adds checking for matching reference and mask in FN6 saves, all while retaining interoperability with FN5 saves.
use std::{path::PathBuf, sync::Mutex};

use crate::sample::ArchivedSample;
use pyo3::prelude::*;
use rayon::prelude::*;

#[pymodule]
/// Fast, efficient SNP distance calculation from disk.
/// Approximately 10x faster than FN5, easier to use and maintain, and adds checking for matching reference and mask in FN6 saves, all while retaining interoperability with FN5 saves.
///
/// This is not intended for comparing large numbers of samples due to performance issues with returning distances rather than writing to stdout. This is intended for computing smaller sets quickly, e.g within an API
fn fn6(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<crate::sample::Sample>()?;
    m.add_class::<crate::sample::SampleHeader>()?;

    m.add_function(wrap_pyfunction!(load_samples, m)?)?;
    m.add_function(wrap_pyfunction!(compute, m)?)?;
    m.add_function(wrap_pyfunction!(crate::sample::distance, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}

/// Load a set of files into memory as Sample objects. This uses multithreading for performance.
///
/// # Arguments
/// - `filepaths`: A vector of paths to load. The file type is determined by the file extension. .fn6 files are expected to be `rkyv` serialized `sample::Sample` structs, while .fn5 files are expected to be in the old format and will be converted to `sample::Sample` structs and then serialized using `rkyv`. FASTA files will be processed and compressed into `sample::Sample` structs and then serialized using `rkyv`.
/// - `reference`: The reference genome sequence as a string. This is only required if at least 1 FASTA file is input
/// - `mask`: A list of positions in the reference genome that should be masked (i.e., ignored) during the analysis. The positions are 0-based. This is only required if at least 1 FASTA file is input.
/// - `mask_hash`: A hash of the mask file. This is used for QC to ensure that the same mask is used for all samples. This is only required if at least 1 FASTA file is input.
/// - `reference_hash`: A hash of the reference genome. This is used for QC to ensure that the same reference is used for all samples. This is only required if at least 1 FASTA file is input.
///
/// # Returns
/// A vector of Sample objects
#[pyfunction(signature = (filepaths, reference = "", mask = Vec::new(), mask_hash = "", reference_hash = ""))]
fn load_samples(
    filepaths: Vec<PathBuf>,
    reference: &str,
    mask: Vec<usize>,
    mask_hash: &str,
    reference_hash: &str,
) -> Vec<crate::sample::Sample> {
    filepaths
        .par_iter()
        .map(|sample_path| {
            let save = crate::load_save(
                sample_path,
                reference,
                &mask,
                mask_hash,
                reference_hash,
            );
            let arch_save = rkyv::access::<ArchivedSample, rkyv::rancor::Error>(&save[..]).unwrap();
            rkyv::deserialize::<crate::sample::Sample, rkyv::rancor::Error>(arch_save).unwrap()
        })
        .collect()
}

/// Given a set of comparisons to do, compute the distances and return them. This uses multithreading for performance, and a mutex to ensure that the output is not interleaved.
///
/// # Arguments
/// - `comparisons`: A vector of tuples, where each tuple contains two byte vectors. Each byte vector is the `rkyv` serialized representation of a `sample::Sample` struct. The order of the byte vectors corresponds to the order of the input filepaths used to load the samples.
/// - `cutoff`: The SNP threshold for distance calculation. If the distance between two samples exceeds this threshold, the distance will not be reported. This is used to speed up distance calculations by allowing for early termination when the distance is large.
fn get_distances(
    comparisons: Vec<(&Vec<u8>, &Vec<u8>)>,
    cutoff: usize,
) -> Vec<(String, String, usize)> {
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
            let dist = crate::sample::arch_distance(sample1, sample2, cutoff);
            if let Some(d) = dist {
                let mut dist_lock = distances.lock().unwrap();
                dist_lock.push((sample1.name.to_string(), sample2.name.to_string(), d));
            }
        })
        .collect::<Vec<()>>();

    // Catch the last bit of distances
    let dist_lock = distances.lock().unwrap();
    dist_lock.clone()
}

/// Compute all distances from a vec of genome save file paths.
///
/// # Arguments
/// - `filepaths`: A vector of files to load. The file type is determined by the file extension. .fn6 files are expected to be `rkyv` serialized `sample::Sample` structs, while .fn5 files are expected to be in the old format and will be converted to `sample::Sample` structs and then serialized using `rkyv`. FASTA files will be reference compressed on the fly and then serialized using `rkyv`.
/// - `reference`: The reference genome sequence as a string. This is only required if at least 1 FASTA file is input
/// - `mask`: A list of positions in the reference genome that should be masked (i.e., ignored) during the analysis. The positions are 0-based. This is only required if at least 1 FASTA file is input.
/// - `mask_hash`: A hash of the mask file. This is used for QC to ensure that the same mask is used for all samples. This is only required if at least 1 FASTA file is input.
/// - `reference_hash`: A hash of the reference genome. This is used for QC to ensure that the same reference is used for all samples. This is only required if at least 1 FASTA file is input.
/// - `cutoff`: The SNP threshold for distance calculation. If the distance between two samples exceeds this threshold, it will not be reported. This is used to speed up distance calculations by allowing for early termination when the distance is large.
/// - `debug`: If true, print debug information about the loading and distance calculation process, including the number of samples loaded, the number of comparisons to be made, and the time taken for each step.
#[pyfunction(signature = (filepaths, reference = "", mask = Vec::new(), mask_hash = "", reference_hash = "", cutoff = 20, debug = false))]
fn compute(
    filepaths: Vec<PathBuf>,
    reference: &str,
    mask: Vec<usize>,
    mask_hash: &str,
    reference_hash: &str,
    cutoff: usize,
    debug: bool,
) -> Vec<(String, String, usize)> {
    // Load the saves
    let start_time = std::time::Instant::now();
    let samples = crate::load_arch_saves(
        filepaths.clone(),
        reference,
        &mask,
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
    let mut distances = Vec::new();
    for (idx, sample1) in samples.iter().enumerate() {
        for sample2 in samples.iter().skip(idx + 1) {
            comparisons.push((sample1, sample2));
            if comparisons.len() > 1_000_000 {
                // We're doing a lot of comparisons, so batch them up to avoid excessive RAM usage
                distances.extend_from_slice(&get_distances(comparisons, cutoff));
                comparisons = Vec::new();
            }
            n_comps += 1;
        }
    }

    // Get last distances
    distances.extend_from_slice(&get_distances(comparisons, cutoff));

    if debug {
        eprintln!(
            "Computed {} distances in {:.2?} ({:.2?}) per comparison",
            n_comps,
            start_time.elapsed(),
            start_time.elapsed() / n_comps as u32
        );
    }

    distances
}
