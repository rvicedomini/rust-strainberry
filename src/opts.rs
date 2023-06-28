use clap::Parser;

#[derive(Parser)]
#[command(version)]
#[command(about = "HiFi-Strainberry: Strain-aware assembly with high-quality long reads", long_about = None)]
pub struct Options {
    
    /// FASTA file of the input assembly
    #[arg(short, long, value_name = "PATH")]
    pub fasta: String,

    /// Long-read alignment in BAM format
    #[arg(short, long, value_name = "PATH")]
    pub bam: String,

    /// Output directory
    #[arg(short, long = "out-dir", value_name = "PATH")]
    pub output_dir: String,
}
