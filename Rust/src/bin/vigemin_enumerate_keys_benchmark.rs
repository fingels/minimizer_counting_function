use clap::Parser;
use std::path::PathBuf;

#[path = "../main.rs"]
#[allow(dead_code)]
mod app;

#[derive(Debug, Parser)]
#[command(name = "vigemin_enumerate_keys_benchmark")]
#[command(about = "Run the optimized Rust benchmark matching scripts/vigemers_enumeration.py")]
struct Cli {
    #[arg(long, default_value_t = 10)]
    m: usize,
    #[arg(long, default_value_t = 31)]
    k: usize,
    #[arg(long)]
    threads: Option<usize>,
    #[arg(long, default_value = "../Figures_theory")]
    output_dir: PathBuf,
    #[arg(long, default_value_t = false)]
    distribution_plots: bool,
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();
    app::run_enumerate_keys(
        cli.m,
        cli.k,
        cli.threads,
        cli.output_dir,
        cli.distribution_plots,
    )
}
