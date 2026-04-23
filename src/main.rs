#![cfg(not(tarpaulin_include))]
use std::path::PathBuf;

use clap::{Parser, Subcommand};
use fn6::sample::{parse_mask, parse_reference};
use rayon::prelude::*;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Args {
    #[command(subcommand)]
    command: Commands,
}
#[derive(Subcommand)]
enum Commands {
    /// Reference compress a sample genome. This will create a .fn6 file that can be used for fast comparisons with other samples. The .fn6 file is a binary file that contains the compressed representation of the sample genome, as well as metadata about the reference and mask used for compression.
    ReferenceCompress {
        /// Path to the reference genome FASTA file
        reference: PathBuf,

        /// Path to the mask file. The mask file is a text file containing the positions of the reference genome that should be masked (i.e., ignored) during the analysis. The positions are 0-based and should be separated by newlines.
        mask: Option<PathBuf>,

        /// Path to the sample genome FASTA file
        sample: PathBuf,

        /// ID for this sample
        #[arg(long)]
        id: Option<String>,

        /// Output path for the .fn6 file. If not provided, the .fn6 file will be saved in the same directory as the sample FASTA file with the same name but with a .fn6 extension.
        #[arg(long)]
        output: Option<PathBuf>,

        /// Whether to print debug information to stderr.
        #[arg(long, default_value_t = false)]
        debug: bool,
    },

    /// Compute SNP distances
    Compute {
        /// Path to the reference genome FASTA file.
        /// Only required if >=1 of the samples specified are fasta files
        #[arg(long, short)]
        reference: Option<PathBuf>,

        /// Path to the mask file. The mask file is a text file containing the positions of the reference genome that should be masked (i.e., ignored) during the analysis. The positions are 0-based and should be separated by newlines.
        /// Only required if >=1 of the samples specified are fasta files
        #[arg(long, short)]
        mask: Option<PathBuf>,

        /// Paths to sample files. Either .fn6, .fn5 or FASTA files.
        /// If FASTA files are provided, the `allow-fasta` flag must also be used. They will then be reference compressed on the fly (using the provided reference and mask) before distance computation.
        #[arg(long, short, num_args = 1..)]
        samples: Option<Vec<PathBuf>>,

        /// Directory to load from. Either .fn6, .fn5 or FASTA files.
        /// If FASTA files are provided, the `allow-fasta` flag must also be used. They will then be reference compressed on the fly (using the provided reference and mask) before distance computation.
        #[arg(long, short)]
        directory: Option<PathBuf>,

        /// SNP threshold
        #[arg(long, default_value_t = 20)]
        cutoff: usize,

        /// Path to a file to store distances in using a space-separated file. Columns are sample1, sample2, distance. If not provided, distances will be printed to stdout in the same format.
        #[arg(long, short)]
        output: Option<PathBuf>,

        /// FASTA file extension to look for when loading from a directory. Only used if loading from a directory and if reference and mask are provided (i.e., if FASTA files need to be reference compressed on the fly). Default is "fasta".
        #[arg(long, default_value = "fasta")]
        fasta_extension: String,

        /// Whether to enable computation from FASTAs. It is recommended to pre-cache the reference compressed versions of the new samples to speed up computation.
        #[arg(long, default_value_t = false)]
        allow_fasta: bool,

        /// Whether to print debug information to stderr.
        #[arg(long, default_value_t = false)]
        debug: bool,
    },

    /// Add some samples to existing samples. Only computes the extra distances required rather than all pairwise distances.
    AddSamples {
        /// Path to the reference genome FASTA file.
        /// Only required if >=1 of the samples specified are fasta files
        #[arg(long, short)]
        reference: Option<PathBuf>,

        /// Path to the mask file. The mask file is a text file containing the positions of the reference genome that should be masked (i.e., ignored) during the analysis. The positions are 0-based and should be separated by newlines.
        /// Only required if >=1 of the samples specified are fasta files
        #[arg(long, short)]
        mask: Option<PathBuf>,

        /// Paths to existing sample files. Either .fn6, .fn5 or FASTA files.
        /// If FASTA files are provided, the `allow-fasta` flag must also be used. They will then be reference compressed on the fly (using the provided reference and mask) before distance computation.
        #[arg(long, short = 's', num_args = 1..)]
        existing_samples: Option<Vec<PathBuf>>,

        /// Directory to load existing saves from. Either .fn6, .fn5 or FASTA files.
        /// If FASTA files are provided, the `allow-fasta` flag must also be used. They will then be reference compressed on the fly (using the provided reference and mask) before distance computation.
        #[arg(long, short = 'd')]
        existing_directory: Option<PathBuf>,

        /// Paths to sample files to add. Either .fn6, .fn5 or FASTA files.
        /// If FASTA files are provided, the `allow-fasta` flag must also be used. They will then be reference compressed on the fly (using the provided reference and mask) before distance computation.
        #[arg(long, short = 'S', num_args = 1..)]
        new_samples: Option<Vec<PathBuf>>,

        /// Directory to load new saves from. Either .fn6, .fn5 or FASTA files.
        /// If FASTA files are provided, the `allow-fasta` flag must also be used. They will then be reference compressed on the fly (using the provided reference and mask) before distance computation.
        #[arg(long, short = 'D')]
        new_directory: Option<PathBuf>,

        /// SNP threshold
        #[arg(long, default_value_t = 20)]
        cutoff: usize,

        /// Path to a file to store distances in using a space-separated file. Columns are sample1, sample2, distance. If not provided, distances will be printed to stdout in the same format.
        #[arg(long, short)]
        output: Option<PathBuf>,

        /// FASTA file extension to look for when loading from a directory. Only used if loading from a directory and if reference and mask are provided (i.e., if FASTA files need to be reference compressed on the fly). Default is "fasta".
        #[arg(long, short, default_value = "fasta")]
        fasta_extension: String,

        /// Whether to enable computation from FASTAs. It is recommended to pre-cache the reference compressed versions of the new samples to speed up computation.
        #[arg(long, default_value_t = false)]
        allow_fasta: bool,

        /// Whether to print debug information to stderr.
        #[arg(long, default_value_t = false)]
        debug: bool,
    },

    /// Reference compress a set of genomes. Dumber than `ReferenceCompress` as it doesn't allow for setting specific IDs or output paths, but much faster as it can be parallelized across samples.
    BulkCompress {
        /// Path to the reference genome FASTA file
        reference: PathBuf,

        /// Path to the mask file. The mask file is a text file containing the positions of the reference genome that should be masked (i.e., ignored) during the analysis. The positions are 0-based and should be separated by newlines.
        mask: Option<PathBuf>,

        /// Paths to sample files. Either FASTA files or .fn6 files
        #[arg(long, short, num_args = 1..)]
        samples: Option<Vec<PathBuf>>,

        /// Directory to load from.
        #[arg(long, short)]
        directory: Option<PathBuf>,

        /// Line separated file to read paths from.
        #[arg(long, short)]
        list: Option<PathBuf>,

        /// Output directory to write saves to. Useful when using `list` as it consolidates saves in a single directory. If not provided, the .fn6 files will be saved in the same directory as their corresponding FASTA files with the same name but with a .fn6 extension.
        #[arg(long, short)]
        output: Option<PathBuf>,

        /// FASTA file extension to look for when loading from a directory
        #[arg(long, default_value = "fasta")]
        fasta_extension: String,

        /// Whether to print debug information to stderr.
        #[arg(long, default_value_t = false)]
        debug: bool,
    },
}

fn main() {
    let args = Args::parse();

    match args.command {
        Commands::ReferenceCompress {
            reference,
            mask,
            sample,
            id,
            output,
            debug,
        } => {
            let reference = parse_reference(reference.as_ref());
            let mask = match mask {
                Some(m) => parse_mask(&m),
                None => Vec::new(),
            };
            if debug {
                eprintln!("Got reference of length {}", reference.len());
                eprintln!("Got mask of length {}", mask.len());
            }

            let reference_hash = sha256::digest(reference.as_bytes());
            let mask_hash = sha256::digest(
                mask.iter()
                    .map(|x| x.to_string())
                    .collect::<String>()
                    .as_bytes(),
            );
            let s = fn6::reference_compress(
                &sample,
                &reference,
                &mask,
                &mask_hash,
                &reference_hash,
                id,
                output,
            );
            // Print the name of the sample and whether it passed QC or not.
            // This enables it to be picked up as/when required by other tools in a pipeline.
            if s.is_qc_passed {
                println!("{}", s.name);
            } else {
                println!("||QC_FAIL: {}||", s.name);
            }
        }
        Commands::BulkCompress {
            reference,
            mask,
            samples,
            directory,
            list,
            output,
            fasta_extension,
            debug,
        } => {
            let start = std::time::Instant::now();
            let reference = parse_reference(reference.as_ref());
            let mask = match mask {
                Some(m) => parse_mask(&m),
                None => Vec::new(),
            };

            let reference_hash = sha256::digest(reference.as_bytes());
            let mask_hash = sha256::digest(
                mask.iter()
                    .map(|x| x.to_string())
                    .collect::<String>()
                    .as_bytes(),
            );

            let mut sample_paths = Vec::new();
            if let Some(samples) = samples {
                sample_paths = samples;
            }
            if let Some(dir) = directory {
                for entry in std::fs::read_dir(dir).unwrap() {
                    let entry = entry.unwrap();
                    let path = entry.path();
                    if path.extension().and_then(|s| s.to_str()) == Some(&fasta_extension) {
                        sample_paths.push(path);
                    }
                }
            }
            if let Some(list) = list {
                let content = std::fs::read_to_string(list).unwrap();
                for line in content.lines() {
                    let path = PathBuf::from(line);
                    if path.extension().and_then(|s| s.to_str()) == Some(&fasta_extension) {
                        sample_paths.push(path);
                    }
                }
            }
            if let Some(output) = output.clone() {
                std::fs::create_dir_all(&output).unwrap();
            }
            let _ = sample_paths
                .par_iter()
                .map(|sample| {
                    let output_path = output
                        .as_ref()
                        .map(|dir| dir.join(sample.file_name().unwrap()).with_extension("fn6"));
                    fn6::reference_compress(
                        sample,
                        &reference,
                        &mask,
                        &mask_hash,
                        &reference_hash,
                        None,
                        output_path,
                    );
                })
                .collect::<Vec<_>>();
            let duration = start.elapsed();
            if debug {
                eprintln!(
                    "Bulk compression completed in {:.2?} seconds",
                    duration.as_secs_f64()
                );
            }
        }
        Commands::Compute {
            samples,
            directory,
            cutoff,
            reference,
            mask,
            fasta_extension,
            allow_fasta,
            debug,
            output,
        } => {
            let mut sample_paths = Vec::new();
            let mut contains_fasta = false;
            if let Some(samples) = samples {
                sample_paths.extend(samples);
            }
            if let Some(dir) = directory {
                for entry in std::fs::read_dir(dir).unwrap() {
                    let entry = entry.unwrap();
                    let path = entry.path();
                    if path.extension().and_then(|s| s.to_str()) == Some("fn6")
                        || path.extension().and_then(|s| s.to_str()) == Some("fn5")
                    {
                        sample_paths.push(path);
                    } else if allow_fasta
                        && path.extension().and_then(|s| s.to_str()) == Some(&fasta_extension)
                    {
                        contains_fasta = true;
                        sample_paths.push(path);
                    }
                }
            }

            if sample_paths.is_empty() {
                eprintln!("No samples provided for distance computation");
                return;
            }

            if !allow_fasta {
                // Fastas aren't allowed so double check we haven't picked any up
                sample_paths.retain(|path| {
                    path.extension().and_then(|s| s.to_str()) == Some("fn6")
                        || path.extension().and_then(|s| s.to_str()) == Some("fn5")
                });
            }

            if contains_fasta && (reference.is_none()) {
                panic!("Reference is required when providing FASTA files");
            }
            let (reference, mask, reference_hash, mask_hash) = match (reference, mask) {
                (Some(r), m) => {
                    let r = parse_reference(r.as_ref());
                    let m = match m {
                        Some(m) => parse_mask(&m),
                        None => Vec::new(),
                    };
                    let reference_hash = sha256::digest(r.as_bytes());
                    let mask_hash = sha256::digest(
                        m.iter()
                            .map(|x| x.to_string())
                            .collect::<String>()
                            .as_bytes(),
                    );
                    (r, m, reference_hash, mask_hash)
                }
                (None, None) => (String::new(), Vec::new(), String::new(), String::new()),
                _ => panic!("Both reference and mask must be provided together"),
            };

            fn6::compute(
                sample_paths,
                &reference,
                &mask,
                &mask_hash,
                &reference_hash,
                cutoff,
                output,
                debug,
            );
        }
        Commands::AddSamples {
            existing_samples,
            existing_directory,
            new_samples,
            new_directory,
            cutoff,
            reference,
            mask,
            fasta_extension,
            allow_fasta,
            debug,
            output,
        } => {
            let mut contains_fasta = false;
            let mut existing_sample_paths = Vec::new();
            if let Some(samples) = existing_samples {
                for sample in samples.iter() {
                    if sample.extension().and_then(|s| s.to_str()) == Some("fn6")
                        || sample.extension().and_then(|s| s.to_str()) == Some("fn5")
                    {
                        existing_sample_paths.push(sample.clone());
                    } else if sample.extension().and_then(|s| s.to_str()) == Some(&fasta_extension)
                    {
                        if allow_fasta {
                            contains_fasta = true;
                            existing_sample_paths.push(sample.clone());
                        } else {
                            eprintln!(
                                "FASTA file passed as input without `--allow-fasta` flag, skipping {}",
                                sample.display()
                            );
                        }
                    }
                }
            }
            if let Some(dir) = existing_directory {
                for entry in std::fs::read_dir(dir).unwrap() {
                    let entry = entry.unwrap();
                    let path = entry.path();
                    if path.extension().and_then(|s| s.to_str()) == Some("fn6")
                        || path.extension().and_then(|s| s.to_str()) == Some("fn5")
                    {
                        existing_sample_paths.push(path);
                    } else if allow_fasta
                        && path.extension().and_then(|s| s.to_str()) == Some(&fasta_extension)
                    {
                        contains_fasta = true;
                        existing_sample_paths.push(path);
                    }
                }
            }

            let mut new_sample_paths = Vec::new();
            if let Some(samples) = new_samples {
                for sample in samples.iter() {
                    if sample.extension().and_then(|s| s.to_str()) == Some("fn6")
                        || sample.extension().and_then(|s| s.to_str()) == Some("fn5")
                    {
                        new_sample_paths.push(sample.clone());
                    } else if sample.extension().and_then(|s| s.to_str()) == Some(&fasta_extension)
                    {
                        if allow_fasta {
                            contains_fasta = true;
                            new_sample_paths.push(sample.clone());
                        } else {
                            eprintln!(
                                "FASTA file passed as input without `--allow-fasta` flag, skipping {}",
                                sample.display()
                            );
                        }
                    }
                }
            }
            if let Some(dir) = new_directory {
                for entry in std::fs::read_dir(dir).unwrap() {
                    let entry = entry.unwrap();
                    let path = entry.path();
                    if path.extension().and_then(|s| s.to_str()) == Some("fn6")
                        || path.extension().and_then(|s| s.to_str()) == Some("fn5")
                    {
                        new_sample_paths.push(path);
                    } else if allow_fasta
                        && path.extension().and_then(|s| s.to_str()) == Some(&fasta_extension)
                    {
                        contains_fasta = true;
                        new_sample_paths.push(path);
                    }
                }
            }

            if new_sample_paths.is_empty() {
                eprintln!("No new samples to add");
                return;
            }

            if contains_fasta && (reference.is_none()) {
                panic!("Reference is required when providing FASTA files");
            }
            let (reference, mask, reference_hash, mask_hash) = match (reference, mask) {
                (Some(r), m) => {
                    let r = parse_reference(r.as_ref());
                    let m = match m {
                        Some(m) => parse_mask(&m),
                        None => Vec::new(),
                    };
                    let reference_hash = sha256::digest(r.as_bytes());
                    let mask_hash = sha256::digest(
                        m.iter()
                            .map(|x| x.to_string())
                            .collect::<String>()
                            .as_bytes(),
                    );
                    (r, m, reference_hash, mask_hash)
                }
                (None, None) => (String::new(), Vec::new(), String::new(), String::new()),
                _ => panic!("Both reference and mask must be provided together"),
            };

            fn6::add_samples(
                existing_sample_paths,
                new_sample_paths,
                &reference,
                &mask,
                &mask_hash,
                &reference_hash,
                cutoff,
                output,
                debug,
            );
        }
    }
}
