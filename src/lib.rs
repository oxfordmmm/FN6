use std::{
    path::{Path, PathBuf},
    sync::Mutex,
};

use crate::sample::ArchivedSample;
use rayon::prelude::*;

pub mod sample;

/// Seemlessly load either a new fasta or an existing save
pub fn load_save(
    filepath: &Path,
    reference: &str,
    mask: &[usize],
    mask_hash: &str,
    reference_hash: &str,
) -> sample::Sample {
    if filepath.to_str().unwrap().to_owned().ends_with(".fn5") {
        let bytes = std::fs::read(filepath).unwrap();
        let arch_sample = unsafe { rkyv::access_unchecked::<ArchivedSample>(&bytes[..]) };
        return rkyv::deserialize::<sample::Sample, rkyv::rancor::Error>(arch_sample).unwrap();
    }
    sample::Sample::new(filepath, reference, mask, mask_hash, reference_hash)
}

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
            path.set_extension("fn5");
            path
        }
    };
    if s.is_qc_passed {
        let serialized = rkyv::to_bytes::<rkyv::rancor::Error>(&s).unwrap();
        std::fs::write(output, serialized).unwrap();
    }

    s
}

pub fn compute(filepaths: Vec<PathBuf>, cutoff: usize) {
    // Load the saves
    let samples: Vec<Vec<u8>> = filepaths
        .par_iter()
        .map(|sample_path| std::fs::read(sample_path).unwrap())
        .collect();

    // Figure out what comparisons we need to do
    let mut comparisons: Vec<(&Vec<u8>, &Vec<u8>)> = Vec::new();
    for (idx, sample1) in samples.iter().enumerate() {
        for sample2 in samples.iter().skip(idx + 1) {
            comparisons.push((sample1, sample2));
        }
    }

    let distances: Mutex<Vec<(String, String, usize)>> = Mutex::new(Vec::new());

    let _ = comparisons
        .par_iter()
        .map(|(sample1, sample2)| {
            let sample1 = unsafe { rkyv::access_unchecked::<ArchivedSample>(&sample1[..]) };
            let sample2 = unsafe { rkyv::access_unchecked::<ArchivedSample>(&sample2[..]) };
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
