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

```bash
cargo run --release -- enumerate-keys --m 10 --k 31 --output-dir ../Figures_theory
```

Optional controls:

- `--threads N` to pin Rayon threads
- `--distribution-plots` to also generate full distribution plots (high memory use):
  - `vigemin_enumeration_k=<k>_m=<m>.png`
  - `vigemin_sorted_enumeration_k=<k>_m=<m>.png`
