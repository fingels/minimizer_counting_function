pub mod dna;
pub mod vigemin;

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use crate::dna::{
    decode_index_to_kmer_inplace, increment_kmer_inplace, parse_dna_word, total_kmers,
};
use crate::vigemin::{KmerCountScratch, VigeminCountingFunction};
use rayon::prelude::*;

const ENUM_CHUNK_SIZE: usize = 1usize << 14;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EnumerationStats {
    pub total_minimizers: u64,
    pub non_zero: u64,
    pub sum_counts: u128,
}

pub fn enumerate_vigemin_counts_parallel(
    m: usize,
    key: &str,
    k: usize,
) -> Result<Vec<u128>, String> {
    let key_codes = parse_dna_word(key)?;
    if key_codes.len() != m {
        return Err(format!(
            "key length must match m, got key length {} and m {}",
            key_codes.len(),
            m
        ));
    }
    if k < m {
        return Err(format!("k must be >= m, got k={k}, m={m}"));
    }

    let total = total_kmers(m)?;
    let total_usize = usize::try_from(total)
        .map_err(|_| format!("4^m={total} does not fit in usize on this platform"))?;
    let mut counts = vec![0u128; total_usize];
    let chunk_size = ENUM_CHUNK_SIZE;

    counts.par_chunks_mut(chunk_size).enumerate().for_each_init(
        || {
            let minimizer_codes = vec![0u8; m];
            let counter =
                VigeminCountingFunction::from_codes_unchecked(&minimizer_codes, &key_codes);
            let scratch = KmerCountScratch::default();
            (minimizer_codes, counter, scratch)
        },
        |(minimizer_codes, counter, scratch), (chunk_idx, chunk)| {
            let start_idx = (chunk_idx * chunk_size) as u64;
            decode_index_to_kmer_inplace(start_idx, minimizer_codes);
            for (offset, slot) in chunk.iter_mut().enumerate() {
                if offset != 0 {
                    increment_kmer_inplace(minimizer_codes);
                }
                counter.reset_from_codes_unchecked(minimizer_codes, &key_codes);
                *slot = counter.kmer_count_unchecked_with_scratch(k, scratch);
            }
        },
    );

    Ok(counts)
}

pub fn enumerate_vigemin_stats_parallel(
    m: usize,
    key: &str,
    k: usize,
) -> Result<EnumerationStats, String> {
    let key_codes = parse_dna_word(key)?;
    if key_codes.len() != m {
        return Err(format!(
            "key length must match m, got key length {} and m {}",
            key_codes.len(),
            m
        ));
    }
    if k < m {
        return Err(format!("k must be >= m, got k={k}, m={m}"));
    }

    let total = total_kmers(m)?;
    let total_usize = usize::try_from(total)
        .map_err(|_| format!("4^m={total} does not fit in usize on this platform"))?;
    let chunk_size = ENUM_CHUNK_SIZE;
    let chunk_count = total_usize.div_ceil(chunk_size);

    let (sum_counts, non_zero) = (0..chunk_count)
        .into_par_iter()
        .map_init(
            || {
                let minimizer_codes = vec![0u8; m];
                let counter =
                    VigeminCountingFunction::from_codes_unchecked(&minimizer_codes, &key_codes);
                let scratch = KmerCountScratch::default();
                (minimizer_codes, counter, scratch)
            },
            |(minimizer_codes, counter, scratch), chunk_idx| {
                let start_idx = chunk_idx * chunk_size;
                let chunk_len = (total_usize - start_idx).min(chunk_size);
                decode_index_to_kmer_inplace(start_idx as u64, minimizer_codes);

                let mut local_sum = 0u128;
                let mut local_non_zero = 0u64;
                for offset in 0..chunk_len {
                    if offset != 0 {
                        increment_kmer_inplace(minimizer_codes);
                    }
                    counter.reset_from_codes_unchecked(minimizer_codes, &key_codes);
                    let count = counter.kmer_count_unchecked_with_scratch(k, scratch);
                    local_sum += count;
                    local_non_zero += u64::from(count > 0);
                }
                (local_sum, local_non_zero)
            },
        )
        .reduce(
            || (0u128, 0u64),
            |(sum_a, nz_a), (sum_b, nz_b)| (sum_a + sum_b, nz_a + nz_b),
        );

    Ok(EnumerationStats {
        total_minimizers: total,
        non_zero,
        sum_counts,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::{decode_index_to_kmer, format_dna_word};

    #[test]
    fn test_known_single_count() {
        let counter = VigeminCountingFunction::new("ACACAA", "CTGGGT").expect("valid inputs");
        let count = counter.kmer_count(10).expect("k valid");
        assert_eq!(count, 31);
    }

    #[test]
    fn test_known_large_count() {
        let counter = VigeminCountingFunction::new("ACACAA", "CTGGGT").expect("valid inputs");
        let count = counter.kmer_count(31).expect("k valid");
        assert_eq!(count, 67_108_863);
    }

    #[test]
    fn test_kmer_counts_up_to_vector() {
        let counter = VigeminCountingFunction::new("AAC", "CAT").expect("valid inputs");
        let values = counter.kmer_counts_up_to(5).expect("k valid");
        assert_eq!(values, vec![1, 7, 21]);
    }

    #[test]
    fn test_enumeration_sum_and_prefix_values() {
        let m = 3;
        let k = 5;
        let key = "CAT";

        let counts = enumerate_vigemin_counts_parallel(m, key, k).expect("valid enumeration");

        let total: u128 = counts.iter().sum();
        assert_eq!(total, 1u128 << (2 * k));

        let expected_first_ten = [9u128, 21, 37, 37, 6, 6, 6, 6, 18, 14];
        assert_eq!(&counts[..10], &expected_first_ten);

        // Spot-check lexicographic minimizer decoding order.
        assert_eq!(format_dna_word(&decode_index_to_kmer(0, 3)), "AAA");
        assert_eq!(format_dna_word(&decode_index_to_kmer(1, 3)), "AAC");
        assert_eq!(format_dna_word(&decode_index_to_kmer(2, 3)), "AAG");
    }

    #[test]
    fn test_enumeration_stats_match_counts() {
        let m = 3;
        let k = 5;
        let key = "CAT";

        let counts = enumerate_vigemin_counts_parallel(m, key, k).expect("valid enumeration");
        let stats = enumerate_vigemin_stats_parallel(m, key, k).expect("valid stats");

        assert_eq!(stats.total_minimizers, counts.len() as u64);
        assert_eq!(stats.sum_counts, counts.iter().sum::<u128>());
        assert_eq!(
            stats.non_zero,
            counts.iter().filter(|&&value| value > 0).count() as u64
        );
    }
}
