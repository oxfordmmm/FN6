# FN6
Fast, efficient and scalable SNP distance calculation from disk.

FN5 reworked into Rust. Approximately 10x faster, easier to use and maintain, and adds checking for matching reference and mask in FN6 saves, all while retaining interoperability with FN5 saves.

# TODO
* Tests


## Usage
All FASTA inputs can be optionally gzipped.
```bash
Usage: fn6 <COMMAND>

Commands:
  reference-compress  Reference compress a sample genome. This will create a .fn6 file that can be used for fast comparisons with other samples. The .fn6 file is a binary file that contains the compressed representation of the sample genome, as well as metadata about the reference and mask used for compression
  compute             Compute distances
  add-samples         Add some samples to existing samples. Only computes the extra distances required rather than all pairwise distances
  bulk-compress       Reference compress a set of genomes. Dumber than `ReferenceCompress` as it doesn't allow for setting specific IDs or output paths, but much faster as it can be parallelized across samples
  help                Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

### Reference compress
```bash
Reference compress a sample genome. This will create a .fn6 file that can be used for fast comparisons with other samples. The .fn6 file is a binary file that contains the compressed representation of the sample genome, as well as metadata about the reference and mask used for compression

Usage: fn6 reference-compress [OPTIONS] <REFERENCE> <MASK> <SAMPLE>

Arguments:
  <REFERENCE>  Path to the reference genome FASTA file
  <MASK>       Path to the mask file. The mask file is a text file containing the positions of the reference genome that should be masked (i.e., ignored) during the analysis. The positions are 0-based and should be separated by newlines
  <SAMPLE>     Path to the sample genome FASTA file

Options:
      --id <ID>          ID for this sample
      --output <OUTPUT>  Output path for the .fn6 file. If not provided, the .fn6 file will be saved in the same directory as the sample FASTA file with the same name but with a .fn6 extension
      --debug            Whether to print debug information to stderr
  -h, --help             Print help

```

### Compute
```bash
Compute SNP distances

Usage: fn6 compute [OPTIONS]

Options:
  -r, --reference <REFERENCE>
          Path to the reference genome FASTA file. Only required if >=1 of the samples specified are fasta files
  -m, --mask <MASK>
          Path to the mask file. The mask file is a text file containing the positions of the reference genome that should be masked (i.e., ignored) during the analysis. The positions are 0-based and should be separated by newlines. Only required if >=1 of the samples specified are fasta files
  -s, --samples <SAMPLES>...
          Paths to sample files. Either .fn6, .fn5 or FASTA files. If FASTA files are provided, the `allow-fasta` flag must also be used. They will then be reference compressed on the fly (using the provided reference and mask) before distance computation
  -d, --directory <DIRECTORY>
          Directory to load from. Either .fn6, .fn5 or FASTA files. If FASTA files are provided, the `allow-fasta` flag must also be used. They will then be reference compressed on the fly (using the provided reference and mask) before distance computation
      --cutoff <CUTOFF>
          SNP threshold [default: 20]
      --fasta-extension <FASTA_EXTENSION>
          FASTA file extension to look for when loading from a directory. Only used if loading from a directory and if reference and mask are provided (i.e., if FASTA files need to be reference compressed on the fly). Default is "fasta" [default: fasta]
      --allow-fasta
          Whether to enable computation from FASTAs. It is recommended to pre-cache the reference compressed versions of the new samples to speed up computation
      --debug
          Whether to print debug information to stderr
  -h, --help
          Print help
```

### Add samples
```bash
Add some samples to existing samples. Only computes the extra distances required rather than all pairwise distances

Usage: fn6 add-samples [OPTIONS]

Options:
  -r, --reference <REFERENCE>
          Path to the reference genome FASTA file. Only required if >=1 of the samples specified are fasta files
  -m, --mask <MASK>
          Path to the mask file. The mask file is a text file containing the positions of the reference genome that should be masked (i.e., ignored) during the analysis. The positions are 0-based and should be separated by newlines. Only required if >=1 of the samples specified are fasta files
  -s, --existing-samples <EXISTING_SAMPLES>...
          Paths to existing sample files. Either .fn6, .fn5 or FASTA files. If FASTA files are provided, the `allow-fasta` flag must also be used. They will then be reference compressed on the fly (using the provided reference and mask) before distance computation
  -d, --existing-directory <EXISTING_DIRECTORY>
          Directory to load existing saves from. Either .fn6, .fn5 or FASTA files. If FASTA files are provided, the `allow-fasta` flag must also be used. They will then be reference compressed on the fly (using the provided reference and mask) before distance computation
  -S, --new-samples <NEW_SAMPLES>...
          Paths to sample files to add. Either .fn6, .fn5 or FASTA files. If FASTA files are provided, the `allow-fasta` flag must also be used. They will then be reference compressed on the fly (using the provided reference and mask) before distance computation
  -D, --new-directory <NEW_DIRECTORY>
          Directory to load new saves from. Either .fn6, .fn5 or FASTA files. If FASTA files are provided, the `allow-fasta` flag must also be used. They will then be reference compressed on the fly (using the provided reference and mask) before distance computation
      --cutoff <CUTOFF>
          SNP threshold [default: 20]
  -f, --fasta-extension <FASTA_EXTENSION>
          FASTA file extension to look for when loading from a directory. Only used if loading from a directory and if reference and mask are provided (i.e., if FASTA files need to be reference compressed on the fly). Default is "fasta" [default: fasta]
      --allow-fasta
          Whether to enable computation from FASTAs. It is recommended to pre-cache the reference compressed versions of the new samples to speed up computation
      --debug
          Whether to print debug information to stderr
  -h, --help
          Print help
```

### Bulk Compress
```bash
Reference compress a set of genomes. Dumber than `ReferenceCompress` as it doesn't allow for setting specific IDs or output paths, but much faster as it can be parallelized across samples

Usage: fn6 bulk-compress [OPTIONS] <REFERENCE> <MASK>

Arguments:
  <REFERENCE>  Path to the reference genome FASTA file
  <MASK>       Path to the mask file. The mask file is a text file containing the positions of the reference genome that should be masked (i.e., ignored) during the analysis. The positions are 0-based and should be separated by newlines

Options:
  -s, --samples <SAMPLES>...
          Paths to sample files. Either FASTA files or .fn6 files
  -d, --directory <DIRECTORY>
          Directory to load from
  -l, --list <LIST>
          Line separated file to read paths from
  -o, --output <OUTPUT>
          Output directory to write saves to. Useful when using `list` as it consolidates saves in a single directory. If not provided, the .fn6 files will be saved in the same directory as their corresponding FASTA files with the same name but with a .fn6 extension
      --fasta-extension <FASTA_EXTENSION>
          FASTA file extension to look for when loading from a directory [default: fasta]
      --debug
          Whether to print debug information to stderr
  -h, --help
          Print help
```

## Performance
Both FN5 and FN6 produce the same results, but with condsiderable time differences. Below is a table of comparison between FN5 and FN6 performance on varying sets of _Mycobacterium tuberculosis_ samples, randomly selected from the CRyPTIC dataset https://doi.org/10.5281/zenodo.16041005
All of the benchmarks were run on the same laptop with an Intel i9-13900H and 32GB RAM, directly from SSD.

| N Samples (passing QC)  | Comparisons | FN5 Reference compression | FN5 Compute pairwise matrix (per comparison) | FN5 Compute pairwise matrix (per comparison) no cutoff | FN6 Reference compression | FN6 Compute pairwise matrix (per comparison) | FN6 Compute pairwise matrix (per comparison) no cutoff | FN6 Compute pairwise matrix from FN5 saves (per comparison) | FN6 Compute pairwise matrix from FN5 saves (per comparison) no cutoff |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| 100         | 4,950     | 1.394s   | 44ms (8.9µs)     | 175ms (35.6µs)     | 0.17s  | 7.93ms (1.60µs)     | 81.28ms (16.42µs) | 9.29ms (1.88µs)     | 53.49ms (10.81µs) |
| 250         | 31,125    | 1.698s   | 208ms (6.7µs)    | 1087ms (34.9µs)    | 0.33s  | 32.29ms (1.04µs)    | 280.36ms (9.01µs) | 24.75ms (795.00ns)  | 270.21ms (8.68µs) |
| 500         | 124,750   | 3.506s   | 628ms (5.0µs)    | 3905ms (31.3µs)    | 0.67s  | 84.04ms (673.00ns)  | 921.34ms (7.39µs) | 48.75ms (390.00ns)  | 1010ms (8.13µs)   |
| 750         | 280,875   | 5.091s   | 1298ms (4.6µs)   | 8179ms (29.1µs)    | 0.92s  | 114.78ms (408.00ns) | 2000ms (7.14µs)   | 91.94ms (327.00ns)  | 5580ms (19.85µs)  |
| 1000        | 499,500   | 6.735s   | 2287ms (4.6µs)   | 14299ms (28.6µs)   | 0.90s  | 152.88ms (306.00ns) | 3480ms (6.98µs)   | 143.89ms (288.00ns) | 8910ms (17.85µs)  |
| 1500 (1493) | 1,113,785 | 11.363s  | 7481ms (6.7µs)   | 100481ms (90.2µs)  | 2.03s  | 300.10ms (269.00ns) | 9130ms (8.20µs)   | 296.32ms (266.00ns) | 22720ms (20.40µs) |
| 2000 (1992) | 1,983,044 | 16.377s  | 30805ms (15.5µs) | 205160ms (103.5µs) | 2.69s  | 487.53ms (245.00ns) | 18390ms (9.27µs)  | 613.04ms (309.00ns) | 42150ms (21.25µs) |

All of the times to compute a pairwise matrix include the time taken to read the saves from disk. The time taken to load the saves from disk varies by tool, but is an important component of the real world runtime.

Mean time to load a save from disk across all above runs:
| FN5 | FN6 | FN6 (from FN5 saves) |
| --- | --- | -------------------- |
| 64.6µs | 25.0µs | 41.3µs |


This also demonstrates the importance of utilising a SNP cutoff within the computation. Otherwise, per sample comparison times scale according to how disparate samples are. 
