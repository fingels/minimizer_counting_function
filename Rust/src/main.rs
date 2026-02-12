use clap::{Parser, Subcommand};
use minimizer_vigemin_rust::dna::{
    decode_index_to_kmer, decode_index_to_kmer_inplace, format_dna_word, parse_dna_word,
    total_kmers,
};
use minimizer_vigemin_rust::vigemin::{KmerCountScratch, VigeminCountingFunction};
use minimizer_vigemin_rust::{enumerate_vigemin_counts_parallel, enumerate_vigemin_stats_parallel};
use plotters::prelude::*;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Debug, Parser)]
#[command(name = "vigemin-fast")]
#[command(about = "High-performance vigemin counting and full minimizer enumeration", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Count how many k-mers admit one specific minimizer for a given key.
    Count {
        #[arg(long)]
        minimizer: String,
        #[arg(long)]
        key: Option<String>,
        #[arg(long)]
        k: usize,
        #[arg(long, default_value_t = false)]
        all_k: bool,
    },
    /// Enumerate counts for all 4^m minimizers for a given key and k.
    Enumerate {
        #[arg(long)]
        m: usize,
        #[arg(long)]
        key: Option<String>,
        #[arg(long)]
        k: usize,
        #[arg(long)]
        threads: Option<usize>,
        #[arg(long)]
        output: Option<PathBuf>,
        #[arg(long, default_value_t = false)]
        sorted: bool,
        #[arg(long, default_value_t = false)]
        non_zero_only: bool,
    },
    /// Reproduce the Python vigemers_enumeration workflow over multiple keys and generate plots.
    EnumerateKeys {
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
    },
    /// Reproduce scripts/pikmin_benchmark.py in Rust, with optional FASTA input.
    Benchmark {
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
    },
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Count {
            minimizer,
            key,
            k,
            all_k,
        } => run_count(&minimizer, key.as_deref(), k, all_k),
        Commands::Enumerate {
            m,
            key,
            k,
            threads,
            output,
            sorted,
            non_zero_only,
        } => run_enumerate(m, key.as_deref(), k, threads, output, sorted, non_zero_only),
        Commands::EnumerateKeys {
            m,
            k,
            threads,
            output_dir,
            distribution_plots,
        } => run_enumerate_keys(m, k, threads, output_dir, distribution_plots),
        Commands::Benchmark {
            k,
            m,
            n_keys,
            seq_size,
            fasta,
            threads,
            seed,
        } => run_benchmark(k, m, n_keys, seq_size, fasta, threads, seed),
    }
}

fn random_dna_key(m: usize) -> String {
    let mut rng = rand::thread_rng();
    random_dna_key_with_rng(m, &mut rng)
}

fn resolve_thread_count(threads: Option<usize>) -> usize {
    threads.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    })
}

const MAX_DENSE_ORACLE_ENTRIES: usize = 1usize << 20;
const MAX_DENSE_COUNTER_ENTRIES: usize = 1usize << 22;

#[derive(Debug, Clone)]
struct PositionBitSet {
    bits: Vec<u64>,
    unique_count: usize,
}

impl PositionBitSet {
    fn new(size: usize) -> Self {
        let words = size.div_ceil(64);
        Self {
            bits: vec![0u64; words],
            unique_count: 0,
        }
    }

    #[inline]
    fn insert(&mut self, position: usize) {
        let word_idx = position >> 6;
        let mask = 1u64 << (position & 63);
        let slot = &mut self.bits[word_idx];
        if (*slot & mask) == 0 {
            *slot |= mask;
            self.unique_count += 1;
        }
    }

    #[inline]
    fn len(&self) -> usize {
        self.unique_count
    }

    #[inline]
    fn contains(&self, position: usize) -> bool {
        let word_idx = position >> 6;
        let mask = 1u64 << (position & 63);
        (self.bits[word_idx] & mask) != 0
    }
}

#[derive(Debug, Clone)]
enum MinimizerCounter {
    Dense(Vec<u64>),
    Sparse(HashMap<u64, u64>),
}

impl MinimizerCounter {
    fn new(total_entries: Option<usize>, dense_threshold: usize) -> Self {
        match total_entries {
            Some(total) if total <= dense_threshold => Self::Dense(vec![0u64; total]),
            _ => Self::Sparse(HashMap::new()),
        }
    }

    #[inline]
    fn increment(&mut self, code: u64) {
        match self {
            Self::Dense(values) => {
                values[code as usize] += 1;
            }
            Self::Sparse(values) => {
                *values.entry(code).or_insert(0) += 1;
            }
        }
    }

    fn non_zero_count(&self) -> usize {
        match self {
            Self::Dense(values) => values.iter().filter(|&&value| value > 0).count(),
            Self::Sparse(values) => values.len(),
        }
    }
}

#[derive(Clone)]
struct LazyOracleLookup {
    key_codes: Vec<u8>,
    minimizer_codes: Vec<u8>,
    counter: VigeminCountingFunction,
    scratch: KmerCountScratch,
    cache: HashMap<u64, u128>,
}

#[derive(Clone)]
enum OracleLookup {
    Dense(Vec<u128>),
    Lazy(Box<LazyOracleLookup>),
}

impl OracleLookup {
    fn new(
        key: &str,
        key_codes: &[u8],
        m: usize,
        k: usize,
        dense_entries: Option<usize>,
    ) -> Result<Self, String> {
        if dense_entries.is_some() {
            let values = enumerate_vigemin_counts_parallel(m, key, k)?;
            return Ok(Self::Dense(values));
        }

        let minimizer_codes = vec![0u8; m];
        let counter = VigeminCountingFunction::from_codes_unchecked(&minimizer_codes, key_codes);
        Ok(Self::Lazy(Box::new(LazyOracleLookup {
            key_codes: key_codes.to_vec(),
            minimizer_codes,
            counter,
            scratch: KmerCountScratch::default(),
            cache: HashMap::new(),
        })))
    }

    #[inline]
    fn get(&mut self, code: u64, k: usize) -> u128 {
        match self {
            Self::Dense(values) => values[code as usize],
            Self::Lazy(state) => {
                if let Some(value) = state.cache.get(&code) {
                    return *value;
                }
                decode_index_to_kmer_inplace(code, &mut state.minimizer_codes);
                state
                    .counter
                    .reset_from_codes_unchecked(&state.minimizer_codes, &state.key_codes);
                let value = state
                    .counter
                    .kmer_count_unchecked_with_scratch(k, &mut state.scratch);
                state.cache.insert(code, value);
                value
            }
        }
    }
}

#[inline]
fn encode_dna_slice(slice: &[u8]) -> u64 {
    let mut value = 0u64;
    for &symbol in slice {
        value = (value << 2) | (symbol as u64);
    }
    value
}

#[inline]
fn find_vigemin_index(kmer: &[u8], key: &[u8]) -> usize {
    let m = key.len();
    let mut best_start = 0usize;
    let max_start = kmer.len() - m;
    for start in 1..=max_start {
        for offset in 0..m {
            let candidate = kmer[start + offset] ^ key[offset];
            let best = kmer[best_start + offset] ^ key[offset];
            if candidate < best {
                best_start = start;
                break;
            }
            if candidate > best {
                break;
            }
        }
    }
    best_start
}

#[inline]
fn find_vigemin_index_and_code(kmer: &[u8], key: &[u8], m: usize) -> (usize, u64) {
    let start = find_vigemin_index(kmer, key);
    let code = encode_dna_slice(&kmer[start..(start + m)]);
    (start, code)
}

fn random_dna_key_with_rng<R: Rng + ?Sized>(m: usize, rng: &mut R) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut key = String::with_capacity(m);
    for _ in 0..m {
        let idx = rng.gen_range(0..bases.len());
        key.push(bases[idx] as char);
    }
    key
}

fn generate_unique_random_keys<R: Rng + ?Sized>(
    m: usize,
    n_keys: usize,
    rng: &mut R,
) -> Result<Vec<String>, String> {
    let total_u64 = total_kmers(m)?;
    let total = usize::try_from(total_u64)
        .map_err(|_| format!("4^m={total_u64} does not fit in usize on this platform"))?;
    if n_keys > total {
        return Err(format!(
            "n_keys={n_keys} is too large for m={m}; maximum distinct keys is {total}"
        ));
    }

    let mut seen = HashSet::with_capacity(n_keys.saturating_mul(2));
    let mut keys = Vec::with_capacity(n_keys);
    while keys.len() < n_keys {
        let idx = rng.gen_range(0..(total as u64));
        if seen.insert(idx) {
            keys.push(format_dna_word(&decode_index_to_kmer(idx, m)));
        }
    }
    Ok(keys)
}

fn generate_random_sequence_codes<R: Rng + ?Sized>(seq_size: usize, rng: &mut R) -> Vec<u8> {
    let mut sequence = Vec::with_capacity(seq_size);
    for _ in 0..seq_size {
        sequence.push(rng.gen_range(0..4u8));
    }
    sequence
}

fn load_fasta_sequence_codes(path: &Path) -> Result<Vec<u8>, String> {
    let file =
        File::open(path).map_err(|e| format!("failed to open FASTA {}: {e}", path.display()))?;
    let reader = BufReader::new(file);

    let mut sequence: Vec<u8> = Vec::new();
    for (line_idx, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|e| {
            format!(
                "failed to read FASTA {} at line {}: {e}",
                path.display(),
                line_idx + 1
            )
        })?;
        if line.starts_with('>') {
            continue;
        }
        for (col_idx, byte) in line.bytes().enumerate() {
            match byte.to_ascii_uppercase() {
                b'A' => sequence.push(0),
                b'C' => sequence.push(1),
                b'G' => sequence.push(2),
                b'T' => sequence.push(3),
                b' ' | b'\t' | b'\r' => {}
                other => {
                    return Err(format!(
                        "invalid FASTA base '{}' at {}:{}:{} (allowed: A/C/G/T)",
                        other as char,
                        path.display(),
                        line_idx + 1,
                        col_idx + 1
                    ));
                }
            }
        }
    }

    if sequence.is_empty() {
        return Err(format!(
            "FASTA {} does not contain any A/C/G/T sequence data",
            path.display()
        ));
    }
    Ok(sequence)
}

#[derive(Clone)]
struct MinQueue {
    positions: Vec<usize>,
    values: Vec<u64>,
    head: usize,
    tail: usize,
    cap: usize,
}

impl MinQueue {
    fn with_window(window_span: usize) -> Self {
        let cap = window_span + 1;
        Self {
            positions: vec![0usize; cap],
            values: vec![0u64; cap],
            head: 0,
            tail: 0,
            cap,
        }
    }

    #[inline]
    fn push(&mut self, pos: usize, value: u64) {
        while self.head != self.tail {
            let last = (self.tail + self.cap - 1) % self.cap;
            if self.values[last] > value {
                self.tail = last;
            } else {
                break;
            }
        }
        self.positions[self.tail] = pos;
        self.values[self.tail] = value;
        self.tail = (self.tail + 1) % self.cap;
    }

    #[inline]
    fn pop_expired(&mut self, min_pos: usize) {
        while self.head != self.tail && self.positions[self.head] < min_pos {
            self.head = (self.head + 1) % self.cap;
        }
    }

    #[inline]
    fn front(&self) -> (usize, u64) {
        (self.positions[self.head], self.values[self.head])
    }
}

struct DenseChunkResult {
    idx: usize,
    key_counts: Vec<u64>,
    heuristic_counts: Vec<u64>,
    heuristic_duplicated_counts: Vec<u64>,
    key_unique_counts: Vec<u64>,
    key_low_masks: Vec<Vec<u64>>,
    key_high_masks: Vec<Vec<u64>>,
    heuristic_unique_count: u64,
    heuristic_low_mask: Vec<u64>,
    heuristic_high_mask: Vec<u64>,
}

#[inline]
fn choose_chunk_windows(window_count: usize, thread_count: usize) -> usize {
    let target_chunks = thread_count.saturating_mul(8).max(1);
    let base = window_count.div_ceil(target_chunks);
    base.max(500_000)
}

#[inline]
fn extract_range_mask(bits: &PositionBitSet, start: usize, len: usize) -> Vec<u64> {
    if len == 0 {
        return Vec::new();
    }
    let mut out = vec![0u64; len.div_ceil(64)];
    for offset in 0..len {
        if bits.contains(start + offset) {
            out[offset >> 6] |= 1u64 << (offset & 63);
        }
    }
    out
}

#[inline]
fn mask_intersection_popcount(left: &[u64], right: &[u64]) -> u64 {
    left.iter()
        .zip(right.iter())
        .map(|(a, b)| u64::from((a & b).count_ones()))
        .sum()
}

#[inline]
fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9E37_79B9_7F4A_7C15);
    value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    value ^ (value >> 31)
}

#[inline]
fn random_base_at(seed: u64, index: usize) -> u8 {
    (splitmix64(seed.wrapping_add(index as u64)) & 0b11) as u8
}

fn process_dense_chunk_from_sequence(
    idx: usize,
    start_window: usize,
    end_window: usize,
    sequence: &[u8],
    k: usize,
    m: usize,
    key_packed_codes: &[u64],
    oracle_dense: &[Vec<u128>],
    total_minimizers: usize,
) -> DenseChunkResult {
    let key_count = key_packed_codes.len();
    let overlap = k - m;
    let window_span = overlap + 1;
    let chunk_windows = end_window - start_window;
    let local_pos_len = chunk_windows + overlap;

    let mut key_counts = vec![0u64; key_count * total_minimizers];
    let mut heuristic_counts = vec![0u64; total_minimizers];
    let mut heuristic_duplicated_counts = vec![0u64; key_count * total_minimizers];

    let mut key_position_sets: Vec<PositionBitSet> = (0..key_count)
        .map(|_| PositionBitSet::new(local_pos_len))
        .collect();
    let mut heuristic_position_set = PositionBitSet::new(local_pos_len);

    let mut queues: Vec<MinQueue> = (0..key_count)
        .map(|_| MinQueue::with_window(window_span))
        .collect();

    let rolling_mask = if 2 * m >= 64 {
        u64::MAX
    } else {
        (1u64 << (2 * m)) - 1
    };

    let mut mmer_code = encode_dna_slice(&sequence[start_window..(start_window + m)]);
    let max_candidate = end_window + overlap - 1;
    for candidate_pos in start_window..=max_candidate {
        if candidate_pos > start_window {
            let next_base = sequence[candidate_pos + m - 1] as u64;
            mmer_code = ((mmer_code << 2) & rolling_mask) | next_base;
        }

        for key_idx in 0..key_count {
            let queue = &mut queues[key_idx];
            let min_valid_pos = candidate_pos.saturating_sub(overlap);
            queue.pop_expired(min_valid_pos);
            let transformed = mmer_code ^ key_packed_codes[key_idx];
            queue.push(candidate_pos, transformed);
        }

        if candidate_pos < start_window + overlap {
            continue;
        }
        let window_start = candidate_pos - overlap;
        if window_start >= end_window {
            break;
        }

        let mut best_oracle = u128::MAX;
        let mut best_key_idx = 0usize;
        let mut best_vigemin_pos = 0usize;
        let mut best_vigemin_code_idx = 0usize;

        for key_idx in 0..key_count {
            let (vigemin_pos, transformed) = queues[key_idx].front();
            let vigemin_code = transformed ^ key_packed_codes[key_idx];
            let vigemin_code_idx = vigemin_code as usize;

            key_counts[key_idx * total_minimizers + vigemin_code_idx] += 1;
            key_position_sets[key_idx].insert(vigemin_pos - start_window);

            let oracle_value = oracle_dense[key_idx][vigemin_code_idx];
            if oracle_value < best_oracle {
                best_oracle = oracle_value;
                best_key_idx = key_idx;
                best_vigemin_pos = vigemin_pos;
                best_vigemin_code_idx = vigemin_code_idx;
            }
        }

        heuristic_counts[best_vigemin_code_idx] += 1;
        heuristic_duplicated_counts[best_key_idx * total_minimizers + best_vigemin_code_idx] += 1;
        heuristic_position_set.insert(best_vigemin_pos - start_window);
    }

    let key_unique_counts = key_position_sets
        .iter()
        .map(|set| set.len() as u64)
        .collect::<Vec<_>>();
    let key_low_masks = key_position_sets
        .iter()
        .map(|set| extract_range_mask(set, 0, overlap))
        .collect::<Vec<_>>();
    let key_high_masks = key_position_sets
        .iter()
        .map(|set| extract_range_mask(set, chunk_windows, overlap))
        .collect::<Vec<_>>();

    DenseChunkResult {
        idx,
        key_counts,
        heuristic_counts,
        heuristic_duplicated_counts,
        key_unique_counts,
        key_low_masks,
        key_high_masks,
        heuristic_unique_count: heuristic_position_set.len() as u64,
        heuristic_low_mask: extract_range_mask(&heuristic_position_set, 0, overlap),
        heuristic_high_mask: extract_range_mask(&heuristic_position_set, chunk_windows, overlap),
    }
}

fn process_dense_chunk_random(
    idx: usize,
    start_window: usize,
    end_window: usize,
    k: usize,
    m: usize,
    key_packed_codes: &[u64],
    oracle_dense: &[Vec<u128>],
    total_minimizers: usize,
    random_seed: u64,
) -> DenseChunkResult {
    let key_count = key_packed_codes.len();
    let overlap = k - m;
    let window_span = overlap + 1;
    let chunk_windows = end_window - start_window;
    let local_pos_len = chunk_windows + overlap;

    let mut key_counts = vec![0u64; key_count * total_minimizers];
    let mut heuristic_counts = vec![0u64; total_minimizers];
    let mut heuristic_duplicated_counts = vec![0u64; key_count * total_minimizers];

    let mut key_position_sets: Vec<PositionBitSet> = (0..key_count)
        .map(|_| PositionBitSet::new(local_pos_len))
        .collect();
    let mut heuristic_position_set = PositionBitSet::new(local_pos_len);

    let mut queues: Vec<MinQueue> = (0..key_count)
        .map(|_| MinQueue::with_window(window_span))
        .collect();

    let rolling_mask = if 2 * m >= 64 {
        u64::MAX
    } else {
        (1u64 << (2 * m)) - 1
    };

    let mut mmer_code = 0u64;
    for pos in start_window..(start_window + m) {
        mmer_code = (mmer_code << 2) | (random_base_at(random_seed, pos) as u64);
    }

    let max_candidate = end_window + overlap - 1;
    if key_count == 2 {
        let key0 = key_packed_codes[0];
        let key1 = key_packed_codes[1];
        let key1_offset = total_minimizers;
        let oracle0 = &oracle_dense[0];
        let oracle1 = &oracle_dense[1];
        let (q0_slice, q1_slice) = queues.split_at_mut(1);
        let q0 = &mut q0_slice[0];
        let q1 = &mut q1_slice[0];
        let (set0_slice, set1_slice) = key_position_sets.split_at_mut(1);
        let set0 = &mut set0_slice[0];
        let set1 = &mut set1_slice[0];

        for candidate_pos in start_window..=max_candidate {
            if candidate_pos > start_window {
                let next_base = random_base_at(random_seed, candidate_pos + m - 1) as u64;
                mmer_code = ((mmer_code << 2) & rolling_mask) | next_base;
            }

            let min_valid_pos = candidate_pos.saturating_sub(overlap);
            q0.pop_expired(min_valid_pos);
            q1.pop_expired(min_valid_pos);
            q0.push(candidate_pos, mmer_code ^ key0);
            q1.push(candidate_pos, mmer_code ^ key1);

            if candidate_pos < start_window + overlap {
                continue;
            }
            let window_start = candidate_pos - overlap;
            if window_start >= end_window {
                break;
            }

            let (vigemin_pos0, transformed0) = q0.front();
            let (vigemin_pos1, transformed1) = q1.front();
            let vigemin_code_idx0 = (transformed0 ^ key0) as usize;
            let vigemin_code_idx1 = (transformed1 ^ key1) as usize;

            key_counts[vigemin_code_idx0] += 1;
            key_counts[key1_offset + vigemin_code_idx1] += 1;
            set0.insert(vigemin_pos0 - start_window);
            set1.insert(vigemin_pos1 - start_window);

            let oracle0_value = oracle0[vigemin_code_idx0];
            let oracle1_value = oracle1[vigemin_code_idx1];
            if oracle1_value < oracle0_value {
                heuristic_counts[vigemin_code_idx1] += 1;
                heuristic_duplicated_counts[key1_offset + vigemin_code_idx1] += 1;
                heuristic_position_set.insert(vigemin_pos1 - start_window);
            } else {
                heuristic_counts[vigemin_code_idx0] += 1;
                heuristic_duplicated_counts[vigemin_code_idx0] += 1;
                heuristic_position_set.insert(vigemin_pos0 - start_window);
            }
        }
    } else {
        for candidate_pos in start_window..=max_candidate {
            if candidate_pos > start_window {
                let next_base = random_base_at(random_seed, candidate_pos + m - 1) as u64;
                mmer_code = ((mmer_code << 2) & rolling_mask) | next_base;
            }

            for key_idx in 0..key_count {
                let queue = &mut queues[key_idx];
                let min_valid_pos = candidate_pos.saturating_sub(overlap);
                queue.pop_expired(min_valid_pos);
                let transformed = mmer_code ^ key_packed_codes[key_idx];
                queue.push(candidate_pos, transformed);
            }

            if candidate_pos < start_window + overlap {
                continue;
            }
            let window_start = candidate_pos - overlap;
            if window_start >= end_window {
                break;
            }

            let mut best_oracle = u128::MAX;
            let mut best_key_idx = 0usize;
            let mut best_vigemin_pos = 0usize;
            let mut best_vigemin_code_idx = 0usize;

            for key_idx in 0..key_count {
                let (vigemin_pos, transformed) = queues[key_idx].front();
                let vigemin_code = transformed ^ key_packed_codes[key_idx];
                let vigemin_code_idx = vigemin_code as usize;

                key_counts[key_idx * total_minimizers + vigemin_code_idx] += 1;
                key_position_sets[key_idx].insert(vigemin_pos - start_window);

                let oracle_value = oracle_dense[key_idx][vigemin_code_idx];
                if oracle_value < best_oracle {
                    best_oracle = oracle_value;
                    best_key_idx = key_idx;
                    best_vigemin_pos = vigemin_pos;
                    best_vigemin_code_idx = vigemin_code_idx;
                }
            }

            heuristic_counts[best_vigemin_code_idx] += 1;
            heuristic_duplicated_counts[best_key_idx * total_minimizers + best_vigemin_code_idx] +=
                1;
            heuristic_position_set.insert(best_vigemin_pos - start_window);
        }
    }

    let key_unique_counts = key_position_sets
        .iter()
        .map(|set| set.len() as u64)
        .collect::<Vec<_>>();
    let key_low_masks = key_position_sets
        .iter()
        .map(|set| extract_range_mask(set, 0, overlap))
        .collect::<Vec<_>>();
    let key_high_masks = key_position_sets
        .iter()
        .map(|set| extract_range_mask(set, chunk_windows, overlap))
        .collect::<Vec<_>>();

    DenseChunkResult {
        idx,
        key_counts,
        heuristic_counts,
        heuristic_duplicated_counts,
        key_unique_counts,
        key_low_masks,
        key_high_masks,
        heuristic_unique_count: heuristic_position_set.len() as u64,
        heuristic_low_mask: extract_range_mask(&heuristic_position_set, 0, overlap),
        heuristic_high_mask: extract_range_mask(&heuristic_position_set, chunk_windows, overlap),
    }
}

fn run_count(minimizer: &str, key: Option<&str>, k: usize, all_k: bool) -> Result<(), String> {
    let start = Instant::now();
    let minimizer_codes = parse_dna_word(minimizer)?;
    let generated_key;
    let key = match key {
        Some(value) => value,
        None => {
            generated_key = random_dna_key(minimizer_codes.len());
            eprintln!("generated_key={generated_key}");
            &generated_key
        }
    };

    let counter = VigeminCountingFunction::new(minimizer, key)?;
    if all_k {
        let values = counter.kmer_counts_up_to(k)?;
        for (offset, value) in values.iter().enumerate() {
            println!("k={}: {}", counter.minimizer_len() + offset, value);
        }
    } else {
        let value = counter.kmer_count(k)?;
        println!("{value}");
    }
    eprintln!("elapsed_ms={}", start.elapsed().as_millis());
    Ok(())
}

pub fn run_benchmark(
    k: usize,
    m: usize,
    n_keys: usize,
    seq_size: usize,
    fasta: Option<PathBuf>,
    threads: Option<usize>,
    seed: Option<u64>,
) -> Result<(), String> {
    if m == 0 {
        return Err("m must be >= 1".to_owned());
    }
    if k < m {
        return Err(format!("k must be >= m, got k={k}, m={m}"));
    }
    if n_keys == 0 {
        return Err("n_keys must be >= 1".to_owned());
    }
    if threads == Some(0) {
        return Err("threads must be >= 1".to_owned());
    }

    let thread_count = resolve_thread_count(threads);
    ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build_global()
        .map_err(|e| format!("failed to set rayon global thread pool: {e}"))?;
    eprintln!("threads={thread_count}");

    let effective_seed = seed.unwrap_or_else(rand::random::<u64>);
    let key_seed = splitmix64(effective_seed ^ 0xC6A4_A793_5BD1_E995);
    let random_sequence_seed = splitmix64(effective_seed ^ 0x9E37_79B9_7F4A_7C15);

    let mut key_rng = StdRng::seed_from_u64(key_seed);
    let keys = generate_unique_random_keys(m, n_keys, &mut key_rng)?;
    println!("k={k} m={m} n_keys={n_keys}");
    println!("seed={effective_seed}");
    println!("keys={}", keys.join(","));

    let key_codes: Vec<Vec<u8>> = keys
        .iter()
        .map(|key| parse_dna_word(key))
        .collect::<Result<Vec<_>, _>>()?;

    let total_minimizers = total_kmers(m)?;
    let total_minimizers_usize = usize::try_from(total_minimizers).ok();
    let dense_oracle_entries = total_minimizers_usize.filter(|&n| n <= MAX_DENSE_ORACLE_ENTRIES);
    let dense_counter_entries = total_minimizers_usize.filter(|&n| n <= MAX_DENSE_COUNTER_ENTRIES);
    let key_packed_codes = key_codes
        .iter()
        .map(|codes| encode_dna_slice(codes))
        .collect::<Vec<_>>();

    let sequence_start = Instant::now();
    let (sequence, sequence_source, sequence_len): (Option<Vec<u8>>, String, usize) =
        if let Some(path) = fasta.as_ref() {
            let loaded = load_fasta_sequence_codes(path)?;
            let length = loaded.len();
            (Some(loaded), format!("fasta:{}", path.display()), length)
        } else if dense_oracle_entries.is_some() {
            (None, format!("random(seq_size={seq_size})"), seq_size)
        } else {
            let mut seq_rng = StdRng::seed_from_u64(random_sequence_seed);
            let generated = generate_random_sequence_codes(seq_size, &mut seq_rng);
            let length = generated.len();
            (
                Some(generated),
                format!("random(seq_size={seq_size})"),
                length,
            )
        };
    let sequence_elapsed = sequence_start.elapsed();
    if sequence_len < k {
        return Err(format!(
            "sequence length must be >= k, got sequence_len={sequence_len}, k={k}"
        ));
    }
    let window_count = sequence_len - k + 1;
    println!("sequence_source={sequence_source}");
    println!("sequence_len={sequence_len} windows={window_count}");
    println!("sequence_build_ms={}", sequence_elapsed.as_millis());

    let oracle_start = Instant::now();
    if let Some(total_minimizers_dense) = dense_oracle_entries {
        let mut oracle_dense = Vec::with_capacity(n_keys);
        for key in &keys {
            oracle_dense.push(enumerate_vigemin_counts_parallel(m, key, k)?);
        }
        let oracle_elapsed = oracle_start.elapsed();
        println!(
            "oracle_mode=dense_precompute oracle_build_ms={}",
            oracle_elapsed.as_millis()
        );

        let chunk_windows = choose_chunk_windows(window_count, thread_count);
        let chunk_count = window_count.div_ceil(chunk_windows);
        println!("scan_chunks={chunk_count} chunk_windows={chunk_windows}");

        let scan_start = Instant::now();
        let mut chunk_results: Vec<DenseChunkResult> = (0..chunk_count)
            .into_par_iter()
            .map(|chunk_idx| {
                let start_window = chunk_idx * chunk_windows;
                let end_window = (start_window + chunk_windows).min(window_count);
                if let Some(sequence_buf) = sequence.as_ref() {
                    process_dense_chunk_from_sequence(
                        chunk_idx,
                        start_window,
                        end_window,
                        sequence_buf,
                        k,
                        m,
                        &key_packed_codes,
                        &oracle_dense,
                        total_minimizers_dense,
                    )
                } else {
                    process_dense_chunk_random(
                        chunk_idx,
                        start_window,
                        end_window,
                        k,
                        m,
                        &key_packed_codes,
                        &oracle_dense,
                        total_minimizers_dense,
                        random_sequence_seed,
                    )
                }
            })
            .collect();
        let scan_elapsed = scan_start.elapsed();

        let merge_start = Instant::now();
        chunk_results.sort_unstable_by_key(|result| result.idx);

        let mut key_partition_counts = vec![0u64; n_keys * total_minimizers_dense];
        let mut heuristic_counts = vec![0u64; total_minimizers_dense];
        let mut heuristic_duplicated_counts = vec![0u64; n_keys * total_minimizers_dense];
        let mut key_unique_totals = vec![0u64; n_keys];
        let mut heuristic_unique_total = 0u64;

        for result in &chunk_results {
            for (dst, src) in key_partition_counts
                .iter_mut()
                .zip(result.key_counts.iter())
            {
                *dst += *src;
            }
            for (dst, src) in heuristic_counts
                .iter_mut()
                .zip(result.heuristic_counts.iter())
            {
                *dst += *src;
            }
            for (dst, src) in heuristic_duplicated_counts
                .iter_mut()
                .zip(result.heuristic_duplicated_counts.iter())
            {
                *dst += *src;
            }
            for (dst, src) in key_unique_totals
                .iter_mut()
                .zip(result.key_unique_counts.iter())
            {
                *dst += *src;
            }
            heuristic_unique_total += result.heuristic_unique_count;
        }

        let overlap = k - m;
        if overlap > 0 {
            for pair in chunk_results.windows(2) {
                let left = &pair[0];
                let right = &pair[1];
                for key_idx in 0..n_keys {
                    key_unique_totals[key_idx] -= mask_intersection_popcount(
                        &left.key_high_masks[key_idx],
                        &right.key_low_masks[key_idx],
                    );
                }
                heuristic_unique_total -= mask_intersection_popcount(
                    &left.heuristic_high_mask,
                    &right.heuristic_low_mask,
                );
            }
        }
        let merge_elapsed = merge_start.elapsed();

        println!(
            "scan_ms={} sec_per_key_kmer={:.9}",
            scan_elapsed.as_millis(),
            scan_elapsed.as_secs_f64() / (window_count as f64 * n_keys as f64)
        );
        println!("merge_ms={}", merge_elapsed.as_millis());

        let random_density = 2.0f64 / ((k - m + 1) as f64);
        let mut key_density_sum = 0.0f64;
        for key_idx in 0..n_keys {
            let start = key_idx * total_minimizers_dense;
            let end = start + total_minimizers_dense;
            let buckets = key_partition_counts[start..end]
                .iter()
                .filter(|&&value| value > 0)
                .count();
            let density = key_unique_totals[key_idx] as f64 / (window_count as f64);
            key_density_sum += density;
            println!(
                "key={} buckets={} density={:.12}",
                keys[key_idx], buckets, density
            );
        }
        println!("Random minimizer density : {random_density}");
        println!(
            "Mean density for the keys : {}",
            key_density_sum / (n_keys as f64)
        );

        let heuristic_bucket_count = heuristic_counts.iter().filter(|&&v| v > 0).count();
        let mut heuristic_duplicated_bucket_count = 0usize;
        for key_idx in 0..n_keys {
            let start = key_idx * total_minimizers_dense;
            let end = start + total_minimizers_dense;
            heuristic_duplicated_bucket_count += heuristic_duplicated_counts[start..end]
                .iter()
                .filter(|&&value| value > 0)
                .count();
        }
        let heuristic_density = heuristic_unique_total as f64 / (window_count as f64);

        println!("Heuristic number of buckets: {heuristic_bucket_count}");
        println!(
            "Heuristic number of buckets (with duplication): {heuristic_duplicated_bucket_count}"
        );
        println!("Random minimizer density : {random_density}");
        println!("Density of the heuristic (duplicated or not): {heuristic_density}");

        println!(
            "elapsed_total_ms={}",
            sequence_elapsed.as_millis()
                + oracle_elapsed.as_millis()
                + scan_elapsed.as_millis()
                + merge_elapsed.as_millis()
        );
        return Ok(());
    }

    let mut oracle_by_key: Vec<OracleLookup> = Vec::with_capacity(n_keys);
    for (key, key_code) in keys.iter().zip(key_codes.iter()) {
        oracle_by_key.push(OracleLookup::new(key, key_code, m, k, None)?);
    }
    let oracle_elapsed = oracle_start.elapsed();
    println!(
        "oracle_mode=lazy_cache oracle_build_ms={}",
        oracle_elapsed.as_millis()
    );

    let mut key_partition_counts: Vec<MinimizerCounter> = (0..n_keys)
        .map(|_| MinimizerCounter::new(dense_counter_entries, MAX_DENSE_COUNTER_ENTRIES))
        .collect();
    let mut key_partition_positions: Vec<PositionBitSet> = (0..n_keys)
        .map(|_| PositionBitSet::new(sequence_len))
        .collect();

    let mut heuristic_counts =
        MinimizerCounter::new(dense_counter_entries, MAX_DENSE_COUNTER_ENTRIES);
    let mut heuristic_duplicated_counts: Vec<MinimizerCounter> = (0..n_keys)
        .map(|_| MinimizerCounter::new(dense_counter_entries, MAX_DENSE_COUNTER_ENTRIES))
        .collect();
    let mut heuristic_positions = PositionBitSet::new(sequence_len);
    let sequence = sequence.as_ref().ok_or_else(|| {
        "internal error: sequence buffer missing for lazy_cache benchmark path".to_owned()
    })?;

    let scan_start = Instant::now();
    for window_start in 0..window_count {
        let kmer = &sequence[window_start..(window_start + k)];

        let mut best_oracle = u128::MAX;
        let mut best_key_idx = 0usize;
        let mut best_vigemin_idx = 0usize;
        let mut best_vigemin_code = 0u64;

        for key_idx in 0..n_keys {
            let (vigemin_idx, vigemin_code) =
                find_vigemin_index_and_code(kmer, &key_codes[key_idx], m);
            key_partition_counts[key_idx].increment(vigemin_code);
            key_partition_positions[key_idx].insert(window_start + vigemin_idx);

            let oracle_value = oracle_by_key[key_idx].get(vigemin_code, k);
            if oracle_value < best_oracle {
                best_oracle = oracle_value;
                best_key_idx = key_idx;
                best_vigemin_idx = vigemin_idx;
                best_vigemin_code = vigemin_code;
            }
        }

        heuristic_counts.increment(best_vigemin_code);
        heuristic_duplicated_counts[best_key_idx].increment(best_vigemin_code);
        heuristic_positions.insert(window_start + best_vigemin_idx);
    }
    let scan_elapsed = scan_start.elapsed();

    println!(
        "scan_ms={} sec_per_key_kmer={:.9}",
        scan_elapsed.as_millis(),
        scan_elapsed.as_secs_f64() / (window_count as f64 * n_keys as f64)
    );

    let random_density = 2.0f64 / ((k - m + 1) as f64);
    let mut key_density_sum = 0.0f64;
    for key_idx in 0..n_keys {
        let density = key_partition_positions[key_idx].len() as f64 / (window_count as f64);
        key_density_sum += density;
        println!(
            "key={} buckets={} density={:.12}",
            keys[key_idx],
            key_partition_counts[key_idx].non_zero_count(),
            density
        );
    }
    println!("Random minimizer density : {random_density}");
    println!(
        "Mean density for the keys : {}",
        key_density_sum / (n_keys as f64)
    );

    let heuristic_bucket_count = heuristic_counts.non_zero_count();
    let heuristic_duplicated_bucket_count: usize = heuristic_duplicated_counts
        .iter()
        .map(MinimizerCounter::non_zero_count)
        .sum();
    let heuristic_density = heuristic_positions.len() as f64 / (window_count as f64);

    println!("Heuristic number of buckets: {heuristic_bucket_count}");
    println!("Heuristic number of buckets (with duplication): {heuristic_duplicated_bucket_count}");
    println!("Random minimizer density : {random_density}");
    println!("Density of the heuristic (duplicated or not): {heuristic_density}");

    println!(
        "elapsed_total_ms={}",
        sequence_elapsed.as_millis() + oracle_elapsed.as_millis() + scan_elapsed.as_millis()
    );
    Ok(())
}

fn run_enumerate(
    m: usize,
    key: Option<&str>,
    k: usize,
    threads: Option<usize>,
    output: Option<PathBuf>,
    sorted: bool,
    non_zero_only: bool,
) -> Result<(), String> {
    if m == 0 {
        return Err("m must be >= 1".to_owned());
    }
    if k < m {
        return Err(format!("k must be >= m, got k={k}, m={m}"));
    }
    if threads == Some(0) {
        return Err("threads must be >= 1".to_owned());
    }

    let generated_key;
    let key = match key {
        Some(value) => value.to_owned(),
        None => {
            generated_key = random_dna_key(m);
            eprintln!("generated_key={generated_key}");
            generated_key
        }
    };

    let key_codes = parse_dna_word(&key)?;
    if key_codes.len() != m {
        return Err(format!(
            "key length must match m, got key length {} and m {}",
            key_codes.len(),
            m
        ));
    }

    let thread_count = resolve_thread_count(threads);
    ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build_global()
        .map_err(|e| format!("failed to set rayon global thread pool: {e}"))?;
    eprintln!("threads={thread_count}");

    let start = Instant::now();
    let (total_minimizers, total_count, non_zero, counts_opt) = if output.is_some() {
        let counts = enumerate_vigemin_counts_parallel(m, &key, k)?;
        let total_minimizers = counts.len() as u64;
        let total_count: u128 = counts.iter().sum();
        let non_zero = counts.iter().filter(|&&v| v > 0).count() as u64;
        (total_minimizers, total_count, non_zero, Some(counts))
    } else {
        let stats = enumerate_vigemin_stats_parallel(m, &key, k)?;
        (
            stats.total_minimizers,
            stats.sum_counts,
            stats.non_zero,
            None,
        )
    };
    let elapsed = start.elapsed();
    let expected_total = (2usize)
        .checked_mul(k)
        .and_then(|shift| 1u128.checked_shl(shift as u32));

    println!("m={m} k={k} key={key}");
    println!("minimizers={total_minimizers}");
    println!("non_zero={non_zero}");
    println!("sum_counts={total_count}");
    if let Some(v) = expected_total {
        println!("expected_total={v}");
    } else {
        println!("expected_total=overflow(u128)");
    }
    println!("elapsed_ms={}", elapsed.as_millis());

    if let Some(path) = output {
        let counts = counts_opt.ok_or_else(|| {
            "internal error: counts missing while output was requested".to_owned()
        })?;
        write_csv(path, &counts, m, sorted, non_zero_only)?;
    }
    Ok(())
}

#[derive(Debug)]
struct KeyEnumeration {
    key: String,
    counts: Vec<u128>,
}

#[derive(Debug)]
struct KeyRunMetrics {
    key: String,
    minimizers: u64,
    non_zero: u64,
    sum_counts: u128,
    elapsed_ms: u128,
    minimizers_per_sec: f64,
}

fn generate_script_keys(m: usize) -> Vec<String> {
    let mut keys = Vec::with_capacity(6);
    keys.push("A".repeat(m));
    keys.push(format!("A{}", "T".repeat(m.saturating_sub(1))));

    let mut alternating = String::with_capacity(m);
    for i in 0..m {
        alternating.push(if i % 2 == 0 { 'A' } else { 'T' });
    }
    keys.push(alternating);

    for prefix in ['C', 'G', 'T'] {
        let mut key = String::with_capacity(m);
        key.push(prefix);
        if m > 1 {
            key.push_str(&random_dna_key(m - 1));
        }
        keys.push(key);
    }

    keys
}

fn log4_count(value: u128) -> f64 {
    if value == 0 {
        -1.0
    } else {
        (value as f64).log(4.0)
    }
}

fn draw_panel<DB: DrawingBackend>(
    area: &DrawingArea<DB, plotters::coord::Shift>,
    records: &[KeyEnumeration],
    color_offset: usize,
    title: &str,
    x_end: i32,
    y_min: f64,
    y_max: f64,
    average_level: f64,
) -> Result<(), String>
where
    DB::ErrorType: std::fmt::Display,
{
    let mut chart = ChartBuilder::on(area)
        .margin(20)
        .caption(title, ("sans-serif", 28))
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(0..(x_end + 1), y_min..y_max)
        .map_err(|e| format!("failed to build chart: {e}"))?;

    chart
        .configure_mesh()
        .x_desc("minimizer index (lex order)")
        .y_desc("log4(count), with zero shown at -1")
        .x_labels(6)
        .y_labels(8)
        .draw()
        .map_err(|e| format!("failed to draw chart mesh: {e}"))?;

    let avg_style = ShapeStyle::from(&RED.mix(0.55)).stroke_width(2);
    let empty_style = ShapeStyle::from(&RED.mix(0.3)).stroke_width(1);
    chart
        .draw_series(LineSeries::new(
            vec![(0, average_level), (x_end, average_level)],
            avg_style,
        ))
        .map_err(|e| format!("failed to draw average line: {e}"))?;
    chart
        .draw_series(LineSeries::new(vec![(0, -1.0), (x_end, -1.0)], empty_style))
        .map_err(|e| format!("failed to draw empty line: {e}"))?;

    for (idx, record) in records.iter().enumerate() {
        let color_idx = color_offset + idx;
        chart
            .draw_series(LineSeries::new(
                record
                    .counts
                    .iter()
                    .enumerate()
                    .map(|(x, &count)| (x as i32, log4_count(count))),
                &Palette99::pick(color_idx),
            ))
            .map_err(|e| format!("failed to draw key series {}: {e}", record.key))?
            .label(format!("gamma={}", record.key))
            .legend(move |(x, y)| {
                PathElement::new(vec![(x, y), (x + 20, y)], Palette99::pick(color_idx))
            });
    }

    chart
        .configure_series_labels()
        .position(SeriesLabelPosition::UpperRight)
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()
        .map_err(|e| format!("failed to draw legend: {e}"))?;

    Ok(())
}

fn plot_throughput(
    metrics: &[KeyRunMetrics],
    m: usize,
    k: usize,
    output_path: &Path,
) -> Result<(), String> {
    if metrics.is_empty() {
        return Err("cannot plot throughput with empty metric set".to_owned());
    }

    let max_throughput = metrics
        .iter()
        .map(|item| item.minimizers_per_sec)
        .fold(0.0_f64, f64::max)
        .max(1.0);

    let root = BitMapBackend::new(output_path, (1800, 1000)).into_drawing_area();
    root.fill(&WHITE)
        .map_err(|e| format!("failed to initialize throughput plot background: {e}"))?;

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .caption(
            format!("Minimizers per second by key (m={m}, k={k})"),
            ("sans-serif", 32),
        )
        .x_label_area_size(120)
        .y_label_area_size(80)
        .build_cartesian_2d(0..(metrics.len() as i32), 0f64..(max_throughput * 1.15))
        .map_err(|e| format!("failed to build throughput chart: {e}"))?;

    chart
        .configure_mesh()
        .x_desc("key index (mapping in CSV/log)")
        .y_desc("minimizers / second")
        .x_labels(metrics.len())
        .y_labels(10)
        .draw()
        .map_err(|e| format!("failed to draw throughput mesh: {e}"))?;

    for (idx, item) in metrics.iter().enumerate() {
        let x0 = idx as i32;
        let x1 = x0 + 1;
        chart
            .draw_series(std::iter::once(Rectangle::new(
                [(x0, 0.0), (x1, item.minimizers_per_sec)],
                Palette99::pick(idx).filled(),
            )))
            .map_err(|e| format!("failed to draw throughput bar for key {}: {e}", item.key))?;
    }

    root.present().map_err(|e| {
        format!(
            "failed to write throughput plot {}: {e}",
            output_path.display()
        )
    })?;
    Ok(())
}

fn write_throughput_csv(path: &Path, metrics: &[KeyRunMetrics]) -> Result<(), String> {
    let file = File::create(path)
        .map_err(|e| format!("failed to create throughput CSV {}: {e}", path.display()))?;
    let mut writer = BufWriter::new(file);
    writer
        .write_all(
            b"key,minimizers,non_zero,sum_counts,elapsed_ms,minimizers_per_second,key_index\n",
        )
        .map_err(|e| format!("failed to write throughput CSV header: {e}"))?;

    for (idx, item) in metrics.iter().enumerate() {
        writeln!(
            writer,
            "{},{},{},{},{},{:.6},{}",
            item.key,
            item.minimizers,
            item.non_zero,
            item.sum_counts,
            item.elapsed_ms,
            item.minimizers_per_sec,
            idx
        )
        .map_err(|e| format!("failed to write throughput CSV row: {e}"))?;
    }

    writer
        .flush()
        .map_err(|e| format!("failed to flush throughput CSV {}: {e}", path.display()))?;
    Ok(())
}

fn plot_unsorted_panels(
    records: &[KeyEnumeration],
    m: usize,
    k: usize,
    output_path: &Path,
) -> Result<(), String> {
    if records.len() < 6 {
        return Err("expected at least 6 key series to plot".to_owned());
    }
    let count_len = records[0].counts.len();
    if count_len == 0 {
        return Err("cannot plot empty series".to_owned());
    }
    let x_end = count_len.saturating_sub(1) as i32;
    let average_level = k.saturating_sub(m) as f64;
    let mut y_max = average_level;
    for record in records {
        for &value in &record.counts {
            y_max = y_max.max(log4_count(value));
        }
    }
    let y_min = -1.2f64;
    let y_max = y_max + 0.8;

    let root = BitMapBackend::new(output_path, (2200, 700)).into_drawing_area();
    root.fill(&WHITE)
        .map_err(|e| format!("failed to initialize plot background: {e}"))?;
    let panels = root.split_evenly((1, 2));

    draw_panel(
        &panels[0],
        &records[..3],
        0,
        "Reference XOR keys",
        x_end,
        y_min,
        y_max,
        average_level,
    )?;
    draw_panel(
        &panels[1],
        &records[3..6],
        3,
        "Randomized XOR keys",
        x_end,
        y_min,
        y_max,
        average_level,
    )?;

    root.present()
        .map_err(|e| format!("failed to write plot {}: {e}", output_path.display()))?;
    Ok(())
}

fn plot_sorted_counts(
    records: &[KeyEnumeration],
    m: usize,
    k: usize,
    output_path: &Path,
) -> Result<(), String> {
    if records.is_empty() {
        return Err("cannot plot empty key list".to_owned());
    }
    let count_len = records[0].counts.len();
    if count_len == 0 {
        return Err("cannot plot empty series".to_owned());
    }
    let x_end = count_len.saturating_sub(1) as i32;
    let average_level = k.saturating_sub(m) as f64;

    let mut y_max = average_level;
    for record in records {
        for &value in &record.counts {
            y_max = y_max.max(log4_count(value));
        }
    }
    let y_min = -1.2f64;
    let y_max = y_max + 0.8;

    let root = BitMapBackend::new(output_path, (1400, 900)).into_drawing_area();
    root.fill(&WHITE)
        .map_err(|e| format!("failed to initialize sorted plot background: {e}"))?;

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .caption(
            format!("Sorted vigemin counts (m={m}, k={k})"),
            ("sans-serif", 30),
        )
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(0..(x_end + 1), y_min..y_max)
        .map_err(|e| format!("failed to build sorted chart: {e}"))?;

    chart
        .configure_mesh()
        .x_desc("sorted minimizer rank")
        .y_desc("log4(count), with zero shown at -1")
        .x_labels(8)
        .y_labels(8)
        .draw()
        .map_err(|e| format!("failed to draw sorted chart mesh: {e}"))?;

    let avg_style = ShapeStyle::from(&RED.mix(0.55)).stroke_width(2);
    let empty_style = ShapeStyle::from(&RED.mix(0.3)).stroke_width(1);
    chart
        .draw_series(LineSeries::new(
            vec![(0, average_level), (x_end, average_level)],
            avg_style,
        ))
        .map_err(|e| format!("failed to draw sorted average line: {e}"))?;
    chart
        .draw_series(LineSeries::new(vec![(0, -1.0), (x_end, -1.0)], empty_style))
        .map_err(|e| format!("failed to draw sorted empty line: {e}"))?;

    for (idx, record) in records.iter().enumerate() {
        let mut sorted = record.counts.clone();
        sorted.sort_unstable_by(|a, b| b.cmp(a));
        let color_idx = idx;
        chart
            .draw_series(LineSeries::new(
                sorted
                    .iter()
                    .enumerate()
                    .map(|(x, &count)| (x as i32, log4_count(count))),
                &Palette99::pick(color_idx),
            ))
            .map_err(|e| format!("failed to draw sorted series {}: {e}", record.key))?
            .label(format!("gamma={}", record.key))
            .legend(move |(x, y)| {
                PathElement::new(vec![(x, y), (x + 20, y)], Palette99::pick(color_idx))
            });
    }

    chart
        .configure_series_labels()
        .position(SeriesLabelPosition::UpperRight)
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()
        .map_err(|e| format!("failed to draw sorted legend: {e}"))?;

    root.present()
        .map_err(|e| format!("failed to write plot {}: {e}", output_path.display()))?;
    Ok(())
}

pub fn run_enumerate_keys(
    m: usize,
    k: usize,
    threads: Option<usize>,
    output_dir: PathBuf,
    distribution_plots: bool,
) -> Result<(), String> {
    if m == 0 {
        return Err("m must be >= 1".to_owned());
    }
    if k < m {
        return Err(format!("k must be >= m, got k={k}, m={m}"));
    }
    if threads == Some(0) {
        return Err("threads must be >= 1".to_owned());
    }

    let thread_count = resolve_thread_count(threads);
    ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build_global()
        .map_err(|e| format!("failed to set rayon global thread pool: {e}"))?;
    eprintln!("threads={thread_count}");

    let keys = generate_script_keys(m);
    println!("m={m} k={k} keys={}", keys.len());
    println!("distribution_plots={distribution_plots}");

    let expected_total = (2usize)
        .checked_mul(k)
        .and_then(|shift| 1u128.checked_shl(shift as u32));
    if expected_total.is_none() {
        eprintln!("expected_total=overflow(u128); sum check disabled");
    }

    let mut metrics: Vec<KeyRunMetrics> = Vec::with_capacity(keys.len());
    let mut distribution_records = if distribution_plots {
        Some(Vec::with_capacity(keys.len()))
    } else {
        None
    };
    let global_start = Instant::now();
    for key in keys {
        let start = Instant::now();
        let (total_minimizers, non_zero, total_count, counts_opt) = if distribution_plots {
            let counts = enumerate_vigemin_counts_parallel(m, &key, k)?;
            let total_minimizers = counts.len() as u64;
            let total_count: u128 = counts.iter().sum();
            let non_zero = counts.iter().filter(|&&v| v > 0).count() as u64;
            (total_minimizers, non_zero, total_count, Some(counts))
        } else {
            let stats = enumerate_vigemin_stats_parallel(m, &key, k)?;
            (
                stats.total_minimizers,
                stats.non_zero,
                stats.sum_counts,
                None,
            )
        };
        let elapsed = start.elapsed();
        let elapsed_ms = elapsed.as_millis();
        let elapsed_s = elapsed.as_secs_f64();
        let minimizers_per_sec = if elapsed_s > 0.0 {
            (total_minimizers as f64) / elapsed_s
        } else {
            total_minimizers as f64
        };

        if let Some(expected) = expected_total {
            if total_count != expected {
                return Err(format!(
                    "sum check failed for key {key}: got {total_count}, expected {expected}"
                ));
            }
        }

        println!(
            "key={key} non_zero={non_zero} sum_counts={total_count} elapsed_ms={elapsed_ms} minimizers_per_sec={minimizers_per_sec:.2}"
        );
        metrics.push(KeyRunMetrics {
            key: key.clone(),
            minimizers: total_minimizers,
            non_zero,
            sum_counts: total_count,
            elapsed_ms,
            minimizers_per_sec,
        });

        if let Some(records) = &mut distribution_records {
            let counts = counts_opt.ok_or_else(|| {
                "internal error: distribution counts missing in distribution mode".to_owned()
            })?;
            records.push(KeyEnumeration { key, counts });
        }
    }
    println!("elapsed_total_ms={}", global_start.elapsed().as_millis());

    fs::create_dir_all(&output_dir).map_err(|e| {
        format!(
            "failed to create output directory {}: {e}",
            output_dir.display()
        )
    })?;

    let throughput_csv = output_dir.join(format!("vigemin_throughput_k={k}_m={m}.csv"));
    let throughput_plot = output_dir.join(format!("vigemin_throughput_k={k}_m={m}.png"));
    write_throughput_csv(&throughput_csv, &metrics)?;
    plot_throughput(&metrics, m, k, &throughput_plot)?;
    eprintln!("wrote {}", throughput_csv.display());
    eprintln!("wrote {}", throughput_plot.display());

    if let Some(records) = distribution_records {
        let unsorted_plot = output_dir.join(format!("vigemin_enumeration_k={k}_m={m}.png"));
        let sorted_plot = output_dir.join(format!("vigemin_sorted_enumeration_k={k}_m={m}.png"));

        plot_unsorted_panels(&records, m, k, &unsorted_plot)?;
        plot_sorted_counts(&records, m, k, &sorted_plot)?;

        eprintln!("wrote {}", unsorted_plot.display());
        eprintln!("wrote {}", sorted_plot.display());
    }

    Ok(())
}

fn write_csv(
    path: PathBuf,
    counts: &[u128],
    m: usize,
    sorted: bool,
    non_zero_only: bool,
) -> Result<(), String> {
    let file = File::create(&path)
        .map_err(|e| format!("failed to create output file {}: {e}", path.display()))?;
    let mut writer = BufWriter::new(file);
    writer
        .write_all(b"minimizer,count\n")
        .map_err(|e| format!("failed to write CSV header: {e}"))?;

    if sorted {
        let mut pairs: Vec<(u64, u128)> = counts
            .iter()
            .enumerate()
            .map(|(idx, &count)| (idx as u64, count))
            .collect();
        pairs.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        for (idx, count) in pairs {
            if non_zero_only && count == 0 {
                continue;
            }
            let minimizer = format_dna_word(&decode_index_to_kmer(idx, m));
            writeln!(writer, "{minimizer},{count}")
                .map_err(|e| format!("failed to write CSV row: {e}"))?;
        }
    } else {
        for (idx, &count) in counts.iter().enumerate() {
            if non_zero_only && count == 0 {
                continue;
            }
            let minimizer = format_dna_word(&decode_index_to_kmer(idx as u64, m));
            writeln!(writer, "{minimizer},{count}")
                .map_err(|e| format!("failed to write CSV row: {e}"))?;
        }
    }

    writer
        .flush()
        .map_err(|e| format!("failed to flush output {}: {e}", path.display()))?;
    eprintln!("wrote {}", path.display());
    Ok(())
}
