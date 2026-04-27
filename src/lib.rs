//! Fast, efficient SNP distance calculation from disk.
//!
//! Approximately 10x faster than FN5, easier to use and maintain, and adds checking for matching reference and mask in FN6 saves, all while retaining interoperability with FN5 saves.
use std::{
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    sync::Mutex,
};

use crate::sample::ArchivedSample;
use rayon::prelude::*;

pub mod py_lib;

pub mod sample;

static MAX_COMPARISONS_IN_MEMORY: usize = 1_000_000;
static MAX_DISTS_IN_MEMORY: usize = 1000;

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
/// - `output`: An optional output path for the distances. If not provided, the distances will be printed to stdout. If provided, the distances will be written to the specified file in the same format as stdout (i.e., "sample1_name sample2_name distance" per line).
///
/// # Output
/// The function prints the distances to stdout as they are computed. Each line of output corresponds to a pairwise comparison and is formatted as "sample1_name sample2_name distance". If the distance exceeds the cutoff, no distance is written.
pub fn get_distances(
    comparisons: Vec<(&Vec<u8>, &Vec<u8>)>,
    cutoff: usize,
    output: Option<PathBuf>,
) {
    let distances: Mutex<Vec<(String, String, usize)>> = Mutex::new(Vec::new());

    let output: Mutex<Box<dyn Write + Send>> = match output {
        Some(path) => Mutex::new(Box::new(BufWriter::new(
            std::fs::File::create(path).unwrap(),
        ))),
        None => Mutex::new(Box::new(BufWriter::new(std::io::stdout()))),
    };

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
                if dist_lock.len() == MAX_DISTS_IN_MEMORY {
                    let mut o = output.lock().unwrap();
                    for (name1, name2, d) in dist_lock.iter() {
                        writeln!(&mut o, "{} {} {}", name1, name2, d).unwrap();
                    }
                    dist_lock.clear();
                }
            }
        })
        .collect::<Vec<()>>();

    // Catch the last bit of distances
    let dist_lock = distances.lock().unwrap();
    let mut o = output.lock().unwrap();
    for (name1, name2, d) in dist_lock.iter() {
        writeln!(&mut o, "{} {} {}", name1, name2, d).unwrap();
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
#[allow(clippy::too_many_arguments)] // It's not that many and they're all important...
pub fn compute(
    filepaths: Vec<PathBuf>,
    reference: &str,
    mask: &[usize],
    mask_hash: &str,
    reference_hash: &str,
    cutoff: usize,
    output: Option<PathBuf>,
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
            if comparisons.len() > MAX_COMPARISONS_IN_MEMORY {
                // We're doing a lot of comparisons, so batch them up to avoid excessive RAM usage
                get_distances(comparisons, cutoff, output.clone());
                comparisons = Vec::new();
            }
            n_comps += 1;
        }
    }

    // Get last distances
    get_distances(comparisons, cutoff, output);

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
    output: Option<PathBuf>,
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
            if comparisons.len() > MAX_COMPARISONS_IN_MEMORY {
                // We're doing a lot of comparisons, so batch them up to avoid excessive RAM usage
                get_distances(comparisons, cutoff, output.clone());
                comparisons = Vec::new();
            }
            n_comps += 1;
        }
    }
    // And each new sample to each other
    for (idx, sample1) in new_samples.iter().enumerate() {
        for sample2 in new_samples.iter().skip(idx + 1) {
            comparisons.push((sample1, sample2));
            if comparisons.len() > MAX_COMPARISONS_IN_MEMORY {
                // We're doing a lot of comparisons, so batch them up to avoid excessive RAM usage
                get_distances(comparisons, cutoff, output.clone());
                comparisons = Vec::new();
            }
            n_comps += 1;
        }
    }

    get_distances(comparisons, cutoff, output);

    if debug {
        eprintln!(
            "Computed {} new distances in {:.2?} ({:.2?}) per comparison",
            n_comps,
            start_time.elapsed(),
            start_time.elapsed() / n_comps as u32
        );
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;

    use super::*;

    macro_rules! assert_panics {
        ($expression:expr) => {
            let result = std::panic::catch_unwind(|| $expression);
            assert!(result.is_err());
        };
    }

    /// Utility function to check if 2 distance files are functionally identical
    /// Useful as there's no guarantee that the pairwise distance (or lines) will be in the same order
    fn compare_dists(expected_dists_path: &str, actual_dists_path: &str) -> bool {
        let expected_output = std::fs::read_to_string(expected_dists_path).unwrap();
        let mut dists1 = expected_output.lines().collect::<Vec<&str>>();
        dists1.sort();

        let actual_output = std::fs::read_to_string(actual_dists_path).unwrap();
        let mut dists2 = actual_output.lines().collect::<Vec<&str>>();
        dists2.sort();

        let dists1 = dists1
            .iter()
            .map(|x| {
                let mut parts = x.split_whitespace().collect::<Vec<&str>>();
                parts.sort();
                parts.join(" ")
            })
            .collect::<Vec<String>>();
        let dists2 = dists2
            .iter()
            .map(|x| {
                let mut parts = x.split_whitespace().collect::<Vec<&str>>();
                parts.sort();
                parts.join(" ")
            })
            .collect::<Vec<String>>();

        dists1 == dists2
    }

    #[test]
    fn test_load_save() {
        // Test loading a non-existent file
        assert_panics!(load_save(
            &PathBuf::from("non_existent_file.fn6"),
            "ACGT",
            &[0],
            "mask_hash",
            "reference_hash"
        ));

        // Test loading a .fn5 file that doesn't exist
        assert_panics!(load_save(
            &PathBuf::from("non_existent_file.fn5"),
            "ACGT",
            &[0],
            "mask_hash",
            "reference_hash"
        ));

        // Test loading a FASTA file that doesn't exist
        assert_panics!(load_save(
            &PathBuf::from("non_existent_file.fasta"),
            "ACGT",
            &[0],
            "mask_hash",
            "reference_hash"
        ));

        // And check real files work as expected
        let reference = sample::parse_reference(Path::new("tests/cases/dummy/reference.fasta"));
        let mask = sample::parse_mask(Path::new("tests/cases/dummy/mask.txt"));
        let reference_hash = sha256::digest(reference.as_bytes());
        let mask_hash = sha256::digest(
            mask.iter()
                .map(|x| x.to_string())
                .collect::<String>()
                .as_bytes(),
        );
        let b_fasta = load_save(
            &PathBuf::from("tests/cases/dummy/1.fasta"),
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
        );
        let b_fn5 = load_save(
            &PathBuf::from("tests/cases/dummy/fn5_saves/uuid1.fn5"),
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
        );
        let b_fn6 = load_save(
            &PathBuf::from("tests/cases/dummy/fn6_saves/1.fn6"),
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
        );
        // Deserialize the saves before comparison for simplicity
        let s_fasta = rkyv::deserialize::<sample::Sample, rkyv::rancor::Error>(
            rkyv::access::<ArchivedSample, rkyv::rancor::Error>(&b_fasta[..]).unwrap(),
        )
        .unwrap();
        let s_fn5 = rkyv::deserialize::<sample::Sample, rkyv::rancor::Error>(
            rkyv::access::<ArchivedSample, rkyv::rancor::Error>(&b_fn5[..]).unwrap(),
        )
        .unwrap();
        let s_fn6 = rkyv::deserialize::<sample::Sample, rkyv::rancor::Error>(
            rkyv::access::<ArchivedSample, rkyv::rancor::Error>(&b_fn6[..]).unwrap(),
        )
        .unwrap();

        // We've implemented PartialEq in such a way that header equality is only checked for
        // cases where it's not an FN5 save (as these don't have headers)
        // So assert_eq should just work here
        assert_eq!(s_fasta, s_fn6);
        assert_eq!(s_fasta, s_fn5);
        assert_eq!(s_fn5, s_fn6);

        // Incorrect or missing references should only cause issues loading from fasta file
        assert_panics!(load_save(
            &PathBuf::from("tests/cases/dummy/1.fasta"),
            "ACGT", // Wrong reference
            &mask,
            &mask_hash,
            &reference_hash
        ));
        assert_panics!(load_save(
            &PathBuf::from("tests/cases/dummy/1.fasta"),
            "", // No reference
            &mask,
            &mask_hash,
            &reference_hash
        ));
        // Not from fn5/6 files, so should panic if reference or mask don't match
        load_save(
            &PathBuf::from("tests/cases/dummy/fn5_saves/uuid1.fn5"),
            "",
            &[],
            "",
            "",
        );
        load_save(
            &PathBuf::from("tests/cases/dummy/fn6_saves/1.fn6"),
            "",
            &[],
            "",
            "",
        );

        // Parsing from fasta with the correct reference but no mask should be fine too
        load_save(
            &PathBuf::from("tests/cases/dummy/1.fasta"),
            &reference,
            &[],
            "",
            &reference_hash,
        );
    }

    #[test]
    fn test_load_arch_saves() {
        // Should be identical to load_save but with multiple files, so we can just check that the outputs are the same for both functions
        let reference = sample::parse_reference(Path::new("tests/cases/dummy/reference.fasta"));
        let mask = sample::parse_mask(Path::new("tests/cases/dummy/mask.txt"));
        let reference_hash = sha256::digest(reference.as_bytes());
        let mask_hash = sha256::digest(
            mask.iter()
                .map(|x| x.to_string())
                .collect::<String>()
                .as_bytes(),
        );
        let filepaths = vec![
            PathBuf::from("tests/cases/dummy/1.fasta"),
            PathBuf::from("tests/cases/dummy/fn5_saves/uuid1.fn5"),
            PathBuf::from("tests/cases/dummy/fn6_saves/1.fn6"),
        ];
        let saves = load_arch_saves(
            filepaths.clone(),
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
        );

        // Deserialize the saves before comparison for simplicity
        let s_fasta = rkyv::deserialize::<sample::Sample, rkyv::rancor::Error>(
            rkyv::access::<ArchivedSample, rkyv::rancor::Error>(&saves[0][..]).unwrap(),
        )
        .unwrap();
        let s_fn5 = rkyv::deserialize::<sample::Sample, rkyv::rancor::Error>(
            rkyv::access::<ArchivedSample, rkyv::rancor::Error>(&saves[1][..]).unwrap(),
        )
        .unwrap();
        let s_fn6 = rkyv::deserialize::<sample::Sample, rkyv::rancor::Error>(
            rkyv::access::<ArchivedSample, rkyv::rancor::Error>(&saves[2][..]).unwrap(),
        )
        .unwrap();

        // We've implemented PartialEq in such a way that header equality is only checked for
        // cases where it's not an FN5 save (as these don't have headers)
        // So assert_eq should just work here
        assert_eq!(s_fasta, s_fn6);
        assert_eq!(s_fasta, s_fn5);
        assert_eq!(s_fn5, s_fn6);
    }

    #[test]
    fn test_reference_compress() {
        let reference = sample::parse_reference(Path::new("tests/cases/dummy/reference.fasta"));
        let mask = sample::parse_mask(Path::new("tests/cases/dummy/mask.txt"));
        let reference_hash = sha256::digest(reference.as_bytes());
        let mask_hash = sha256::digest(
            mask.iter()
                .map(|x| x.to_string())
                .collect::<String>()
                .as_bytes(),
        );

        let output_path = PathBuf::from("tests/output/dummy-1.fn6");
        if output_path.exists() {
            std::fs::remove_file(&output_path).unwrap();
        }
        let s = reference_compress(
            &PathBuf::from("tests/cases/dummy/1.fasta"),
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
            Some("test_sample".to_string()),
            Some(output_path.clone()),
        );
        assert!(s.is_qc_passed);
        assert!(output_path.exists());

        // Sample which fails QC shouldn't get an output file
        let output_path2 = PathBuf::from("tests/cases/dummy/qc-fail.fn6");
        if output_path2.exists() {
            std::fs::remove_file(&output_path2).unwrap();
        }
        let s = reference_compress(
            &PathBuf::from("tests/cases/dummy/qc-fail.fasta"),
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
            None,
            None,
        );
        assert!(!s.is_qc_passed);
        assert!(!output_path2.exists());
    }

    #[test]
    fn test_get_distances() {
        let reference = sample::parse_reference(Path::new("tests/cases/dummy/reference.fasta"));
        let mask = sample::parse_mask(Path::new("tests/cases/dummy/mask.txt"));
        let reference_hash = sha256::digest(reference.as_bytes());
        let mask_hash = sha256::digest(
            mask.iter()
                .map(|x| x.to_string())
                .collect::<String>()
                .as_bytes(),
        );
        let b_fasta = load_save(
            &PathBuf::from("tests/cases/dummy/1.fasta"),
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
        );
        let b_fn5 = load_save(
            &PathBuf::from("tests/cases/dummy/fn5_saves/uuid1.fn5"),
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
        );
        let b_fn6 = load_save(
            &PathBuf::from("tests/cases/dummy/fn6_saves/1.fn6"),
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
        );
        let b2_fn6 = load_save(
            &PathBuf::from("tests/cases/dummy/fn6_saves/2.fn6"),
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
        );

        // We know these samples are identical, so distance should be 0
        get_distances(
            vec![(&b_fasta, &b_fn5), (&b_fasta, &b_fn6), (&b_fn5, &b_fn6)],
            10,
            Some(PathBuf::from("tests/output/dummy_distances.txt")),
        );
        let output = std::fs::read_to_string("tests/output/dummy_distances.txt").unwrap();
        let lines = output.lines().collect::<Vec<&str>>();
        assert_eq!(lines.len(), 0);

        // The 2.fn6 sample has a SNP at position 0, so distance should be 1 to everything else
        get_distances(
            vec![(&b_fasta, &b2_fn6), (&b_fn5, &b2_fn6), (&b_fn6, &b2_fn6)],
            10,
            Some(PathBuf::from("tests/output/dummy_distances2.txt")),
        );
        let output = std::fs::read_to_string("tests/output/dummy_distances2.txt").unwrap();
        let lines = output.lines().collect::<Vec<&str>>();
        assert_eq!(lines.len(), 3);
        for line in lines {
            let parts = line.split_whitespace().collect::<Vec<&str>>();
            assert_eq!(parts.len(), 3);
            let dist: usize = parts[2].parse().unwrap();
            assert_eq!(dist, 1);
        }
    }

    #[test]
    fn test_compute() {
        let reference = sample::parse_reference(Path::new("tests/cases/dummy/reference.fasta"));
        let mask = sample::parse_mask(Path::new("tests/cases/dummy/mask.txt"));
        let reference_hash = sha256::digest(reference.as_bytes());
        let mask_hash = sha256::digest(
            mask.iter()
                .map(|x| x.to_string())
                .collect::<String>()
                .as_bytes(),
        );

        // Just check that it runs without error and produces the expected output file
        let output_path = PathBuf::from("tests/output/dummy_compute_distances1.txt");
        if output_path.exists() {
            std::fs::remove_file(&output_path).unwrap();
        }
        compute(
            vec![
                PathBuf::from("tests/cases/dummy/fn6_saves/1.fn6"),
                PathBuf::from("tests/cases/dummy/fn6_saves/2.fn6"),
                PathBuf::from("tests/cases/dummy/fn6_saves/3.fn6"),
                PathBuf::from("tests/cases/dummy/fn6_saves/4.fn6"),
                PathBuf::from("tests/cases/dummy/fn6_saves/5.fn6"),
            ],
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
            20,
            Some(output_path.clone()),
            false,
        );
        assert!(output_path.exists());
        assert!(compare_dists(
            "tests/cases/dummy/all-distances.txt",
            output_path.to_str().unwrap(),
        ));

        // Should be unchanged if we run it again without reference/mask
        let output_path = PathBuf::from("tests/output/dummy_compute_distances2.txt");
        if output_path.exists() {
            std::fs::remove_file(&output_path).unwrap();
        }
        compute(
            vec![
                PathBuf::from("tests/cases/dummy/fn6_saves/1.fn6"),
                PathBuf::from("tests/cases/dummy/fn6_saves/2.fn6"),
                PathBuf::from("tests/cases/dummy/fn6_saves/3.fn6"),
                PathBuf::from("tests/cases/dummy/fn6_saves/4.fn6"),
                PathBuf::from("tests/cases/dummy/fn6_saves/5.fn6"),
            ],
            "",
            &[],
            "",
            "",
            20,
            Some(output_path.clone()),
            false,
        );
        assert!(output_path.exists());
        assert!(compare_dists(
            "tests/cases/dummy/all-distances.txt",
            output_path.to_str().unwrap(),
        ));

        // FN5 saves should work
        let output_path = PathBuf::from("tests/output/dummy_compute_distances3.txt");
        if output_path.exists() {
            std::fs::remove_file(&output_path).unwrap();
        }
        compute(
            vec![
                PathBuf::from("tests/cases/dummy/fn5_saves/uuid1.fn5"),
                PathBuf::from("tests/cases/dummy/fn5_saves/uuid2.fn5"),
                PathBuf::from("tests/cases/dummy/fn5_saves/uuid3.fn5"),
                PathBuf::from("tests/cases/dummy/fn5_saves/uuid4.fn5"),
                PathBuf::from("tests/cases/dummy/fn5_saves/uuid5.fn5"),
            ],
            "",
            &[],
            "",
            "",
            20,
            Some(output_path.clone()),
            false,
        );
        assert!(output_path.exists());
        assert!(compare_dists(
            "tests/cases/dummy/all-distances.txt",
            output_path.to_str().unwrap(),
        ));

        // As should direct FASTA input
        let output_path = PathBuf::from("tests/output/dummy_compute_distances4.txt");
        if output_path.exists() {
            std::fs::remove_file(&output_path).unwrap();
        }
        compute(
            vec![
                PathBuf::from("tests/cases/dummy/1.fasta"),
                PathBuf::from("tests/cases/dummy/2.fasta"),
                PathBuf::from("tests/cases/dummy/3.fasta"),
                PathBuf::from("tests/cases/dummy/4.fasta"),
                PathBuf::from("tests/cases/dummy/5.fasta"),
            ],
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
            20,
            Some(output_path.clone()),
            false,
        );
        assert!(output_path.exists());
        assert!(compare_dists(
            "tests/cases/dummy/all-distances.txt",
            output_path.to_str().unwrap(),
        ));

        // Adding a sample which fails QC shouldn't change distances as it should not be considered
        let output_path = PathBuf::from("tests/output/dummy_compute_distances5.txt");
        if output_path.exists() {
            std::fs::remove_file(&output_path).unwrap();
        }
        compute(
            vec![
                PathBuf::from("tests/cases/dummy/1.fasta"),
                PathBuf::from("tests/cases/dummy/2.fasta"),
                PathBuf::from("tests/cases/dummy/3.fasta"),
                PathBuf::from("tests/cases/dummy/4.fasta"),
                PathBuf::from("tests/cases/dummy/5.fasta"),
                PathBuf::from("tests/cases/dummy/qc-fail.fasta"),
            ],
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
            20,
            Some(output_path.clone()),
            false,
        );
        assert!(output_path.exists());
        assert!(compare_dists(
            "tests/cases/dummy/all-distances.txt",
            output_path.to_str().unwrap(),
        ));

        // Mixing any of the above should be fine too
        let output_path = PathBuf::from("tests/output/dummy_compute_distances6.txt");
        if output_path.exists() {
            std::fs::remove_file(&output_path).unwrap();
        }
        compute(
            vec![
                PathBuf::from("tests/cases/dummy/1.fasta"),
                PathBuf::from("tests/cases/dummy/fn6_saves/2.fn6"),
                PathBuf::from("tests/cases/dummy/fn6_saves/3.fn6"),
                PathBuf::from("tests/cases/dummy/fn6_saves/4.fn6"),
                PathBuf::from("tests/cases/dummy/fn5_saves/uuid5.fn5"),
                PathBuf::from("tests/cases/dummy/qc-fail.fasta"),
            ],
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
            20,
            Some(output_path.clone()),
            false,
        );
        assert!(output_path.exists());
        assert!(compare_dists(
            "tests/cases/dummy/all-distances.txt",
            output_path.to_str().unwrap(),
        ));

        // Cutoff should work as expected
        // All dists here are 1 except one distance which is 2, so with a cutoff of 1 we should get all but that distance
        let output_path = PathBuf::from("tests/output/dummy_compute_distances6.txt");
        if output_path.exists() {
            std::fs::remove_file(&output_path).unwrap();
        }
        compute(
            vec![
                PathBuf::from("tests/cases/dummy/1.fasta"),
                PathBuf::from("tests/cases/dummy/fn6_saves/2.fn6"),
                PathBuf::from("tests/cases/dummy/fn6_saves/3.fn6"),
                PathBuf::from("tests/cases/dummy/fn6_saves/4.fn6"),
                PathBuf::from("tests/cases/dummy/fn5_saves/uuid5.fn5"),
                PathBuf::from("tests/cases/dummy/qc-fail.fasta"),
            ],
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
            1,
            Some(output_path.clone()),
            false,
        );
        assert!(output_path.exists());
        assert!(compare_dists(
            "tests/cases/dummy/all-distances-cutoff1.txt",
            output_path.to_str().unwrap(),
        ));
        // Similarly, cutoff of 0 shouldn't return any distances
        let output_path = PathBuf::from("tests/output/dummy_compute_distances6.txt");
        if output_path.exists() {
            std::fs::remove_file(&output_path).unwrap();
        }
        compute(
            vec![
                PathBuf::from("tests/cases/dummy/1.fasta"),
                PathBuf::from("tests/cases/dummy/fn6_saves/2.fn6"),
                PathBuf::from("tests/cases/dummy/fn6_saves/3.fn6"),
                PathBuf::from("tests/cases/dummy/fn6_saves/4.fn6"),
                PathBuf::from("tests/cases/dummy/fn5_saves/uuid5.fn5"),
                PathBuf::from("tests/cases/dummy/qc-fail.fasta"),
            ],
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
            0,
            Some(output_path.clone()),
            false,
        );
        assert!(output_path.exists());
        assert!(compare_dists(
            "tests/cases/dummy/no-distances.txt",
            output_path.to_str().unwrap(),
        ));
    }

    #[test]
    fn test_add_samples() {
        let reference = sample::parse_reference(Path::new("tests/cases/dummy/reference.fasta"));
        let mask = sample::parse_mask(Path::new("tests/cases/dummy/mask.txt"));
        let reference_hash = sha256::digest(reference.as_bytes());
        let mask_hash = sha256::digest(
            mask.iter()
                .map(|x| x.to_string())
                .collect::<String>()
                .as_bytes(),
        );

        // Just check that it runs without error and produces the expected output file
        let output_path = PathBuf::from("tests/output/dummy_add_samples_distances1.txt");
        if output_path.exists() {
            std::fs::remove_file(&output_path).unwrap();
        }
        add_samples(
            vec![
                PathBuf::from("tests/cases/dummy/fn6_saves/1.fn6"),
                PathBuf::from("tests/cases/dummy/2.fasta"),
            ],
            vec![
                PathBuf::from("tests/cases/dummy/fn6_saves/3.fn6"),
                PathBuf::from("tests/cases/dummy/fn6_saves/4.fn6"),
                PathBuf::from("tests/cases/dummy/fn5_saves/uuid5.fn5"),
            ],
            &reference,
            &mask,
            &mask_hash,
            &reference_hash,
            20,
            Some(output_path.clone()),
            false,
        );
        assert!(output_path.exists());
        assert!(compare_dists(
            "tests/cases/dummy/adding-samples.txt",
            output_path.to_str().unwrap(),
        ));
    }
}
