use clap::Parser;
use std::path::PathBuf;

#[path = "../main.rs"]
#[allow(dead_code)]
mod app;

#[derive(Debug, Parser)]
#[command(name = "pikmin_benchmark")]
#[command(about = "Run the optimized Rust implementation of scripts/pikmin_benchmark.py")]
struct Cli {
    #[arg(long, default_value_t = 31)]
    k: usize,
    #[arg(long, default_value_t = 6)]
    m: usize,
    #[arg(long, default_value_t = 8)]
    n_keys: usize,
    #[arg(long, default_value_t = 100_000)]
    seq_size: usize,
    #[arg(long)]
    fasta: Option<PathBuf>,
    #[arg(long)]
    threads: Option<usize>,
    #[arg(long)]
    seed: Option<u64>,
    #[arg(long, default_value = ".")]
    output_dir: PathBuf,
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();
    app::run_benchmark(
        cli.k,
        cli.m,
        cli.n_keys,
        cli.seq_size,
        cli.fasta,
        cli.threads,
        cli.seed,
        cli.output_dir,
    )
}
