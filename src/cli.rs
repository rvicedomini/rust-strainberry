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

    /// Maximum number of theads
    #[arg(short = 't', long = "threads", value_name = "NUM", default_value_t = 1)]
    pub nb_threads: usize,

    /// Stop after phasing and read partitioning
    #[arg(long = "phase-only")]
    pub phase_only: bool,

    /// Do not split at putative misassembly events
    #[arg(long = "no-split")]
    pub no_split: bool,

    /// Minimum QUAL value for loaded variants (effective only with --vcf)
    #[arg(long = "min-var-qual", value_name = "NUM", default_value_t = 0)]
    pub min_var_qual: usize,

    /// Minimum number of phased variants to retain a haplotype
    #[arg(long = "min-snv", value_name = "NUM", default_value_t = 3)]
    pub min_snv: usize,

    /// Minimum number of shared positions between a sequence and a haplotype to consider a match
    #[arg(long = "min-shared-snv", value_name = "NUM", default_value_t = 1)]
    pub min_shared_snv: usize,

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
    #[arg(long = "min-overhang", value_name = "NUM", default_value_t = 500)]
    pub min_overhang: usize,

    /// Minimum aware-contig length
    #[arg(long = "min-aware-ctg-len", value_name = "NUM", default_value_t = 2000)]
    pub min_aware_ctg_len: usize,

}


// TODO: validate parsed options and possibly estimate other parameters
pub struct Config;

impl Config {

    pub fn from_options(_opts:Options) -> Config {
        todo!()
    }
}
