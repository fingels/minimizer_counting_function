# Rust Vigemin Engine

High-performance Rust rewrite of the vigemin counting core, with a parallel enumeration binary.

## Build

```bash
cargo build --release
```

Release profile is configured for aggressive optimization:

- `opt-level = 3`
- `lto = "fat"`
- `codegen-units = 1`
- `panic = "abort"`
- `target-cpu = native` via `.cargo/config.toml`

## CLI

### 1. Count one minimizer

```bash
cargo run --release -- count --minimizer ACACAA --key CTGGGT --k 10
```

`--key` is optional. If omitted, a random valid DNA key of length `m` is generated and printed to stderr.

```bash
cargo run --release -- count --minimizer ACACAA --k 10
```

Print all values from `k = m` to a target:

```bash
cargo run --release -- count --minimizer AAC --key CAT --k 5 --all-k
```

### 2. Enumerate all minimizers (`4^m`) in parallel

```bash
cargo run --release -- enumerate --m 10 --key CTGGGT --k 31
```

Optional controls:

- `--threads N` to pin Rayon threads
- `--output path.csv` to write CSV (`minimizer,count`)
- `--sorted` to sort descending by count before CSV write
- `--non-zero-only` to omit empty minimizers in CSV

Defaults:

- If `--key` is omitted, a random valid DNA key of length `m` is generated and printed to stderr.
- If `--threads` is omitted, all available CPU threads are used.

Example:

```bash
cargo run --release -- enumerate --m 3 --key CAT --k 5 --output counts.csv --sorted
```

### 3. Reproduce `vigemers_enumeration.py` workflow (multi-key + plots)

This command loops over the same key family as the Python script:

- `A...A` (lexicographic)
- `A` + `T...T` (anti-lexicographic pattern)
- alternating `ATAT...`
- three random keys prefixed by `C`, `G`, `T`

Then it computes per-key stats, verifies sums, and writes throughput artifacts:

- `vigemin_throughput_k=<k>_m=<m>.csv`
- `vigemin_throughput_k=<k>_m=<m>.png`

Subcommand form:

```bash
cargo run --release -- enumerate-keys --m 10 --k 31 --output-dir ../Figures_theory
```

Dedicated binary form:

```bash
cargo run --release --bin vigemin_enumerate_keys_benchmark -- --m 10 --k 31 --output-dir ../Figures_theory
```

Optional controls:

- `--threads N` to pin Rayon threads
- `--distribution-plots` to also generate full distribution plots (high memory use):
  - `vigemin_enumeration_k=<k>_m=<m>.png`
  - `vigemin_sorted_enumeration_k=<k>_m=<m>.png`

### 4. Reproduce `pikmin_benchmark.py` workflow (optimized scan + heuristic)

This command computes, in one pass:

- per-key vigemin partitions
- heuristic partition (min oracle across keys)
- heuristic partition with duplicated key/minimizer buckets
- densities and bucket counts

Default behavior matches the Python benchmark setup with random DNA sequence generation:

Subcommand form:

```bash
cargo run --release -- benchmark --k 31 --m 6 --n-keys 8 --seq-size 100000
```

Dedicated binary form:

```bash
cargo run --release --bin pikmin_benchmark -- --k 31 --m 6 --n-keys 8 --seq-size 100000
```

Use an input FASTA file instead of a random sequence:

```bash
cargo run --release --bin pikmin_benchmark -- --k 31 --m 6 --n-keys 8 --fasta /path/to/input.fa
```

Optional controls:

- `--threads N` to pin Rayon threads (used by oracle precompute path)
- `--seed S` for deterministic random key generation (and deterministic random sequence if `--fasta` is not used)
  If `--seed` is omitted, a random seed is generated and printed.
- `--output-dir PATH` directory where the benchmark plot is written (default: current directory)

This command writes:

- `pikmin_benchmark_k=<k>_m=<m>_N_keys=<n_keys>_seq_size=<sequence_len>.png`

## Benchmark Binaries

### `vigemin_enumerate_keys_benchmark`

Runs the `vigemers_enumeration.py`-style benchmark (multi-key throughput + optional distribution plots).

```bash
cargo run --release --bin vigemin_enumerate_keys_benchmark -- --m 10 --k 31 --output-dir ../Figures_theory
```

Options:

- `--m <usize>` default: `10`
- `--k <usize>` default: `31`
- `--threads <usize>` default: all available threads
- `--output-dir <path>` default: `../Figures_theory`
- `--distribution-plots` default: off

### `pikmin_benchmark`

Runs the optimized Rust version of `scripts/pikmin_benchmark.py`.

```bash
cargo run --release --bin pikmin_benchmark -- --k 31 --m 6 --n-keys 8 --seq-size 100000
```

FASTA mode:

```bash
cargo run --release --bin pikmin_benchmark -- --k 31 --m 6 --n-keys 8 --fasta /path/to/input.fa
```

Options:

- `--k <usize>` default: `31`
- `--m <usize>` default: `6`
- `--n-keys <usize>` default: `8`
- `--seq-size <usize>` default: `100000` (used only when `--fasta` is not provided)
- `--fasta <path>` optional; if set, input sequence is read from FASTA
- `--threads <usize>` default: all available threads
- `--seed <u64>` optional; if omitted, a random seed is generated and printed
- `--output-dir <path>` default: `.`
