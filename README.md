# Strainberry 2

Strainberry2 is the new and improved version of Strainberry. It is conceived to perform strain-level metagenome assembly with modern long reads (e.g., PacBio HiFi, ONT R10.4).


## Installation

### Prerequisites

- Rust >= 1.82

### Software dependencies

Strainberry2 needs the following two bioinformatics tool to be available:
- samtools (tested with version 1.15)
- minimap2 (tested with version 2.28)

### Build from source

Clone the repository:

```
git clone https://github.com/rvicedomini/rust-strainberry.git
```

To build the executable:

```
cd rust-strainberry && RUSTFLAGS="-C target-cpu=native" cargo build --release
```

The `strainberry2` executable will be available in the directory `./target/release/`


## Quick usage

Strainberry2 requires only two input files:
- A reference assembly, which could be provided in either FASTA or GFA format.
- A set of long reads in FASTQ or FASTA format.

All input files could be provided either uncompressed or gzipped.

To run Strainberry2 with **ONT R10.4** reads using 8 threads:
```
strainberry2 --in-ont reads.fastq.gz -r assembly.fasta --out-dir output -t 8
```

To run Strainberry2 with **PacBio HiFi** reads using 8 threads:
```
strainberry2 --in-hifi reads.fastq.gz -r assembly.fasta --out-dir output -t 8
```


## Advanced usage

Coming soon...

## Output files

The strain-aware output assembly is written in the `assembly.fasta` file in inside the output directory provided with the `--out-dir` parameter. Moreover, an assembly graph (in GFA format) is written in the `assembly.gfa` file.

## Example

Coming soon...


