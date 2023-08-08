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

    /// User provided vcf file of SNV positions to consider
    #[arg(long = "vcf", value_name = "PATH")]
    pub vcf_file: Option<String>,

    /// Lookback distance
    #[arg(short = 'l', long = "lookback", value_name = "NUM", default_value_t = 3000)]
    pub lookback: usize,

    /// Minimum MAPQ value to consider a read alignment
    #[arg(short = 'q', long = "min-mapq", value_name = "NUM", default_value_t = 20)]
    pub min_mapq: u8,

    /// Minimum overhang length
    #[arg(short = 't', long = "threads", value_name = "NUM", default_value_t = 1)]
    pub nb_threads: usize,

    /// Minimum number of alternative-allele observations
    #[arg(long = "min-alt-count", value_name = "NUM", default_value_t = 5)]
    pub min_alt_count: usize,

    /// Minimum fraction of alternative-allele observations
    #[arg(long = "min-alt-frac", value_name = "FLOAT", default_value_t = 0.125)]
    pub min_alt_frac: f64,

    /// Minimum indel length
    #[arg(long = "min-indel", value_name = "NUM", default_value_t = 100)]
    pub min_indel: usize,

    /// Minimum overhang length
    #[arg(long = "min-overhang", value_name = "NUM", default_value_t = 1000)]
    pub min_overhang: usize,

    /// Minimum aware-contig length
    #[arg(long = "min-aware-ctg-len", value_name = "NUM", default_value_t = 3000)]
    pub min_aware_ctg_len: usize,

}


// TODO: validate parsed options and possibly estimate other parameters
pub struct Config;

impl Config {

    pub fn from_options(_opts:Options) -> Config {
        todo!()
    }
}
