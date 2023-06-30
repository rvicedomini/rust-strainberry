use clap::Parser;

#[derive(Parser)]
#[command(version)]
#[command(about = "HiFi-Strainberry: Strain-aware assembly with high-quality long reads", long_about = None)]
pub struct Options {
    
    /// FASTA file of the input assembly
    #[arg(short = 'f', long = "fasta", value_name = "PATH")]
    pub fasta_file: String,

    /// Long-read alignment in BAM format
    #[arg(short = 'b', long = "bam", value_name = "PATH")]
    pub bam_file: String,

    /// Output directory
    #[arg(short = 'o', long = "out-dir", value_name = "PATH")]
    pub output_dir: String,

    /// Minimum number of consecutive k-mers to define a marker
    #[arg(long = "vcf", value_name = "PATH")]
    pub vcf_file: Option<String>,

    /// Lookback distance
    #[arg(short = 'l', long = "lookback", value_name = "NUM")]
    pub lookback: Option<usize>,

    /// Minimum number of consecutive k-mers to define a marker
    #[arg(short = 'q', long = "min-mapq", value_name = "NUM", default_value_t = 20)]
    pub min_mapq: u8,

    /// Minimum number of consecutive k-mers to define a marker
    #[arg(long = "min-alt-count", value_name = "NUM", default_value_t = 5)]
    pub min_alt_count: usize,

    /// Relative difference between the distance of a k-mer on the reference/query to be included in a marker
    #[arg(long = "min-alt-frac", value_name = "FLOAT", default_value_t = 0.125)]
    pub min_alt_frac: f64,
}