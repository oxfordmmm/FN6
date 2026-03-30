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
        mask: PathBuf,

        /// Path to the sample genome FASTA file
        sample: PathBuf,

        /// ID for this sample
        #[arg(long)]
        id: Option<String>,

        /// Output path for the .fn6 file. If not provided, the .fn6 file will be saved in the same directory as the sample FASTA file with the same name but with a .fn6 extension.
        #[arg(long)]
        output: Option<PathBuf>,
    },

    /// Compute distances
    Compute {
        /// Paths to sample files. Either FASTA files or .fn6 files
        #[arg(long, short, num_args = 1..)]
        samples: Option<Vec<PathBuf>>,

        /// Directory to load from.
        #[arg(long, short)]
        directory: Option<PathBuf>,

        /// SNP threshold
        #[arg(long, default_value_t = 20)]
        cutoff: usize,
    },

    /// Reference compress a set of genomes. Dumber than `ReferenceCompress` as it doesn't allow for setting specific IDs or output paths, but much faster as it can be parallelized across samples.
    BulkCompress {
        /// Path to the reference genome FASTA file
        reference: PathBuf,

        /// Path to the mask file. The mask file is a text file containing the positions of the reference genome that should be masked (i.e., ignored) during the analysis. The positions are 0-based and should be separated by newlines.
        mask: PathBuf,

        /// Paths to sample files. Either FASTA files or .fn6 files
        #[arg(long, short, num_args = 1..)]
        samples: Option<Vec<PathBuf>>,

        /// Directory to load from.
        #[arg(long, short)]
        directory: Option<PathBuf>,
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
        } => {
            let reference = parse_reference(reference.as_ref());
            eprintln!("Got reference of length {}", reference.len());
            let mut m = Vec::new();
            let mask = parse_mask(mask.as_ref(), &mut m);
            eprintln!("Got mask of length {}", mask.len());

            let reference_hash = sha256::digest(reference.as_bytes());
            let mask_hash = sha256::digest(
                mask.iter()
                    .map(|x| x.to_string())
                    .collect::<String>()
                    .as_bytes(),
            );
            fn6::reference_compress(
                &sample,
                &reference,
                mask,
                &mask_hash,
                &reference_hash,
                id,
                output,
            );
        }
        Commands::BulkCompress {
            reference,
            mask,
            samples,
            directory,
        } => {
            let reference = parse_reference(reference.as_ref());
            eprintln!("Got reference of length {}", reference.len());
            let mut m = Vec::new();
            let mask = parse_mask(mask.as_ref(), &mut m);
            eprintln!("Got mask of length {}", mask.len());

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
                    if path.extension().and_then(|s| s.to_str()) == Some("fasta")
                        || path.extension().and_then(|s| s.to_str()) == Some("fna")
                        || path.extension().and_then(|s| s.to_str()) == Some("fa")
                    {
                        sample_paths.push(path);
                    }
                }
            }
            eprintln!("Got {} samples to compress", sample_paths.len());
            let _ = sample_paths
                .par_iter()
                .map(|sample| {
                    fn6::reference_compress(
                        sample,
                        &reference,
                        mask,
                        &mask_hash,
                        &reference_hash,
                        None,
                        None,
                    )
                })
                .collect::<Vec<_>>();
        }
        Commands::Compute {
            samples,
            directory,
            cutoff,
        } => {
            let mut sample_paths = Vec::new();
            if let Some(samples) = samples {
                sample_paths = samples
                    .into_iter()
                    .filter(|p| {
                        p.extension().and_then(|s| s.to_str()) == Some("fn6")
                            || p.extension().and_then(|s| s.to_str()) == Some("fn5")
                    })
                    .collect();
            }
            if let Some(dir) = directory {
                for entry in std::fs::read_dir(dir).unwrap() {
                    let entry = entry.unwrap();
                    let path = entry.path();
                    if path.extension().and_then(|s| s.to_str()) == Some("fn6")
                        || path.extension().and_then(|s| s.to_str()) == Some("fn5")
                    {
                        sample_paths.push(path);
                    }
                }
            }
            fn6::compute(sample_paths, cutoff);
        }
    }
}
