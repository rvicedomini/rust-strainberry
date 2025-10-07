# Strainberry 2

Strainberry2 is the new and improved version of Strainberry. It is conceived to perform strain-level metagenome assembly with modern long reads (e.g., PacBio HiFi, ONT R10.4).


## Installation

Strainberry has been developed and tested under a Linux environment.

### Prerequisites

- Rust >= 1.82

### Software dependencies

Strainberry2 requires the following two bioinformatics tool to be available in the system:
- samtools >= 1.15
- minimap2 >= 2.27

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

Strainberry2 requires two input files:
- **A reference assembly**, which could be provided in either FASTA or GFA format.
- **A set of long reads**, possibly in FASTQ format (FASTA is accepted but providing quality values might improve final assembly, especially with ONT data)

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

The strain-aware output assembly is written in the `assembly.fasta` file in inside the output directory provided with the `--out-dir` parameter. A corresponding assembly graph in GFA format is also available in the `assembly.gfa` file.

## Example

In order to verify that Strainberry2 has been correctly installed, it is possible to test it on a tiny dataset in the `example` sub-directory.

Given a strain-oblivious assembly (e.g., file `assembly.fa.gz`) and a collection of accurate ONT reads (e.g., file `reads.fa.gz`), it is possible to run Strainberry2 using 4 threads as follows:

```
cd example
strainberry2 --in-ont reads.fa.gz -r assembly.fa.gz -o output -t 4
```

If everything worked as expected, you will find the strain-aware assembly `assembly.fasta` and the corresponding strain-aware assembly graph `assembly.gfa` in the `output` directory. Both files should contain 5 contigs, one for each of the 5 *Escherichia coli* strain sequences present in the input data. 



