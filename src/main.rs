use clap::Parser;
use fn5::sample::{Sample, distance, parse_mask, parse_reference};
use rayon::prelude::*;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Path to the reference genome FASTA file
    reference: String,

    /// Path to the mask file. The mask file is a text file containing the positions of the reference genome that should be masked (i.e., ignored) during the analysis. The positions are 0-based and should be separated by newlines.
    mask: String,

    /// Paths to the sample genome FASTA files
    #[arg(long, required = true, num_args = 1..)]
    samples: Vec<String>,

    /// SNP threshold
    #[arg(long, default_value_t = 20)]
    cutoff: usize,
}

fn main() {
    let args = Args::parse();
    let reference = parse_reference(args.reference.as_ref());
    eprintln!("Got reference of length {}", reference.len());
    let mut m = Vec::new();
    let mask = parse_mask(args.mask.as_ref(), &mut m);
    eprintln!("Got mask of length {}", mask.len());

    let reference_hash = sha256::digest(args.reference.as_bytes());
    let mask_hash = sha256::digest(args.mask.as_bytes());

    let samples: Vec<Sample> = args
        .samples
        .par_iter()
        .map(|sample_path| {
            Sample::new(
                sample_path.as_ref(),
                reference.as_str(),
                mask,
                &mask_hash,
                &reference_hash,
            )
        })
        .collect();

    eprintln!(
        "Got samples {:?}",
        samples
    );

    let mut comparisons: Vec<(&Sample, &Sample)> = Vec::new();
    for (idx1, sample1) in samples.iter().enumerate() {
        for idx2 in idx1 + 1..samples.len() {
            let sample2 = &samples[idx2];
            comparisons.push((sample1, sample2));
        }
    }
    
    let _ = comparisons.par_iter().map(|(sample1, sample2)| {
        if sample1.name == sample2.name {
            return;
        }
        let dist = distance(sample1, sample2, args.cutoff);
        match dist {
            Some(d) => println!("{} {} {}", sample1.name, sample2.name, d),
            None => return,
        }
    }).collect::<Vec<()>>();
}
