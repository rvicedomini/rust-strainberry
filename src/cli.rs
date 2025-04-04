use clap::{Parser, ValueEnum};

#[derive(Parser)]
#[command(version)]
#[command(about = "HiFi-Strainberry: Strain-aware assembly with high-quality long reads", long_about = None)]
pub struct Options {
    
    /// Input assembly in FASTA format
    #[arg(short = 'f', long = "fasta", value_name = "PATH")]
    pub fasta_file: String,

    /// Input HiFi reads
    #[arg(long = "in-hifi", value_name = "PATH")]
    pub in_hifi: Option<String>,
    
    /// Input ONT reads
    #[arg(long = "in-ont", value_name = "PATH")]
    pub in_ont: Option<String>,

    /// Long-read alignment in BAM format
    #[arg(short = 'b', long = "bam", value_name = "PATH")]
    pub bam_file: String,

    /// User provided vcf file of SNV positions to consider
    #[arg(short = 'v', long = "vcf", value_name = "PATH")]
    pub vcf_file: Option<String>,

    /// Output directory
    #[arg(short = 'o', long = "out-dir", value_name = "PATH")]
    pub output_dir: String,

    /// Lookback distance
    #[arg(short = 'l', long = "lookback", value_name = "NUM", default_value_t = 3000)]
    pub lookback: usize,

    /// Minimum MAPQ value to consider a read alignment
    #[arg(short = 'q', long = "min-mapq", value_name = "NUM", default_value_t = 20)]
    pub min_mapq: u8,

    /// Maximum number of theads
    #[arg(short = 't', long = "threads", value_name = "NUM", default_value_t = 1)]
    pub nb_threads: usize,

    /// Do not delete temporary files
    #[arg(long = "keep-temp")]
    pub keep_temp: bool,

    /// Disable post-assembly polishing
    #[arg(long = "no-polish")]
    pub no_polish: bool,

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

    /// Sequencing read technology
    #[arg(value_enum, short='m', long="mode", value_name="STR", default_value_t = Mode::Hifi)]
    pub mode: Mode,
}


#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Mode {
    Hifi,
    Nano,
}


// TODO: validate parsed options and possibly estimate other parameters
pub struct Config;

impl Config {

    pub fn from_options(_opts:Options) -> Config {
        todo!()
    }
}
