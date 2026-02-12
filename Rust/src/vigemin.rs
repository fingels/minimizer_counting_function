use crate::dna::{
    ALL_DNA_MASK, DNA_ALPHABET_SIZE, dna_mask, dna_mask_len, lex_rank, parse_dna_word, xor_symbol,
};

#[derive(Clone, Copy, Eq, PartialEq)]
enum CmpTag {
    Eq,
    Lt,
    Gt,
}

#[derive(Default, Clone)]
pub struct KmerCountScratch {
    antemer_prefix: Vec<u128>,
    antemer_array: Vec<u128>,
    postmer_prefix: Vec<u128>,
    postmer_array: Vec<u128>,
    postmer_values: Vec<u128>,
}

#[derive(Clone)]
pub struct VigeminCountingFunction {
    length: usize,
    antemer_max_prefix_size: usize,
    postmer_max_size: usize,
    prefix_letters_vectors: Vec<[usize; DNA_ALPHABET_SIZE]>,
    autocorrelation_matrix: Vec<CmpTag>,
    suffix_key_convolution: Vec<bool>,
    antemer_alphabet_zero: Vec<u8>,
    antemer_alphabet_not_zero: Vec<u8>,
    postmer_small_values_alphabet_zero: Vec<u8>,
    postmer_small_values_alphabet_not_zero: Vec<u8>,
    postmer_high_values_alphabet: Vec<Vec<(usize, u8)>>,
    alphabet_i: Vec<u8>,
}

impl VigeminCountingFunction {
    pub fn new(minimizer: &str, key: &str) -> Result<Self, String> {
        let minimizer_codes = parse_dna_word(minimizer)?;
        let key_codes = parse_dna_word(key)?;
        Self::from_codes(&minimizer_codes, &key_codes)
    }

    pub fn from_codes(minimizer: &[u8], key: &[u8]) -> Result<Self, String> {
        if minimizer.is_empty() {
            return Err("minimizer cannot be empty".to_string());
        }
        if minimizer.len() != key.len() {
            return Err(format!(
                "key and minimizer must have the same length, got {} and {}",
                key.len(),
                minimizer.len()
            ));
        }
        Ok(Self::from_codes_unchecked(minimizer, key))
    }

    #[inline]
    pub fn from_codes_unchecked(minimizer: &[u8], key: &[u8]) -> Self {
        debug_assert!(!minimizer.is_empty());
        debug_assert_eq!(minimizer.len(), key.len());

        let length = minimizer.len();
        let inf = usize::MAX / 4;

        let mut antemer_max_prefix_size = length;
        let mut postmer_max_size = usize::MAX;

        let mut prefix_letters_vectors = vec![[inf; DNA_ALPHABET_SIZE]; length + 1];
        let mut autocorrelation_matrix = vec![CmpTag::Eq; length * length];

        let mut alphabet_i = vec![0u8; length + 1];
        for i in 0..length {
            let ref_xor = minimizer[i] ^ key[i];
            let mut mask = 0u8;
            for a in 0..DNA_ALPHABET_SIZE {
                let cand = (a as u8) ^ key[i];
                if cand > ref_xor {
                    mask |= dna_mask(a as u8);
                }
            }
            alphabet_i[i] = mask;
        }
        alphabet_i[length] = ALL_DNA_MASK;

        prefix_letters_vectors[1] = [0; DNA_ALPHABET_SIZE];
        prefix_letters_vectors[1][minimizer[0] as usize] = 2;

        for j in 1..length {
            let mut flag = false;
            let mut symb = CmpTag::Eq;

            for i in j..length {
                if !flag {
                    let t = i - j;
                    let left = minimizer[i] ^ key[t];
                    let right = minimizer[t] ^ key[t];
                    if left < right {
                        symb = CmpTag::Lt;
                        flag = true;
                        antemer_max_prefix_size = antemer_max_prefix_size.min(i + 1);
                        postmer_max_size = postmer_max_size.min(j + 1);
                    } else if left > right {
                        symb = CmpTag::Gt;
                        flag = true;
                    } else {
                        symb = CmpTag::Eq;
                        let letter = minimizer[t + 1] as usize;
                        prefix_letters_vectors[i + 1][letter] =
                            prefix_letters_vectors[i + 1][letter].min(j + 1);
                    }
                }
                autocorrelation_matrix[i * length + j] = symb;
            }

            for a in 0..DNA_ALPHABET_SIZE {
                if prefix_letters_vectors[j + 1][a] == inf {
                    prefix_letters_vectors[j + 1][a] =
                        if (a as u8) == minimizer[0] { j + 2 } else { 0 };
                }
            }
        }

        let mut suffix_key_convolution = vec![false; antemer_max_prefix_size];
        for i in 0..antemer_max_prefix_size {
            suffix_key_convolution[i] = suffix_key_convolution_gt(&minimizer, key, i);
        }

        let mut antemer_alphabet_zero = vec![0u8; length + 1];
        let mut antemer_alphabet_not_zero = vec![0u8; length + 1];

        for i in 1..length {
            let mut correct_alphabet = alphabet_i[i] & (alphabet_i[0] | dna_mask(minimizer[0]));

            for j in 1..i {
                if autocorrelation_matrix[(i - 1) * length + j] == CmpTag::Eq {
                    correct_alphabet &= alphabet_i[i - j] | dna_mask(minimizer[i - j]);
                }
            }

            for a in 0..DNA_ALPHABET_SIZE {
                let bit = dna_mask(a as u8);
                if correct_alphabet & bit == 0 {
                    continue;
                }
                if prefix_letters_vectors[i][a] == 0 {
                    antemer_alphabet_zero[i] |= bit;
                } else {
                    antemer_alphabet_not_zero[i] |= bit;
                }
            }
        }

        let mut postmer_small_values_alphabet_zero = vec![0u8; length + 1];
        let mut postmer_small_values_alphabet_not_zero = vec![0u8; length + 1];
        let mut postmer_high_values_alphabet = vec![Vec::new(); length + 1];

        for i in 1..=length {
            for a in 0..DNA_ALPHABET_SIZE {
                let bit = dna_mask(a as u8);
                if prefix_letters_vectors[i][a] == 0 {
                    postmer_small_values_alphabet_zero[i] |= bit;
                } else {
                    postmer_small_values_alphabet_not_zero[i] |= bit;
                }
            }

            if i < length {
                let minimizer_bit = dna_mask(minimizer[i]);
                postmer_small_values_alphabet_zero[i] &= !minimizer_bit;
                postmer_small_values_alphabet_not_zero[i] &= !minimizer_bit;
            }

            postmer_high_values_alphabet[i].push((i, alphabet_i[0] | dna_mask(minimizer[0])));
            for k in 1..i {
                if autocorrelation_matrix[(i - 1) * length + k] == CmpTag::Eq {
                    postmer_high_values_alphabet[i]
                        .push((k, alphabet_i[i - k] | dna_mask(minimizer[i - k])));
                }
            }
        }

        Self {
            length,
            antemer_max_prefix_size,
            postmer_max_size,
            prefix_letters_vectors,
            autocorrelation_matrix,
            suffix_key_convolution,
            antemer_alphabet_zero,
            antemer_alphabet_not_zero,
            postmer_small_values_alphabet_zero,
            postmer_small_values_alphabet_not_zero,
            postmer_high_values_alphabet,
            alphabet_i,
        }
    }

    #[inline]
    pub fn reset_from_codes_unchecked(&mut self, minimizer: &[u8], key: &[u8]) {
        debug_assert!(!minimizer.is_empty());
        debug_assert_eq!(minimizer.len(), key.len());

        if self.length != minimizer.len() {
            *self = Self::from_codes_unchecked(minimizer, key);
            return;
        }

        let length = self.length;
        let inf = usize::MAX / 4;

        self.antemer_max_prefix_size = length;
        self.postmer_max_size = usize::MAX;
        self.prefix_letters_vectors.fill([inf; DNA_ALPHABET_SIZE]);
        self.autocorrelation_matrix.fill(CmpTag::Eq);
        self.antemer_alphabet_zero.fill(0);
        self.antemer_alphabet_not_zero.fill(0);
        self.postmer_small_values_alphabet_zero.fill(0);
        self.postmer_small_values_alphabet_not_zero.fill(0);
        self.alphabet_i.fill(0);
        for values in &mut self.postmer_high_values_alphabet {
            values.clear();
        }

        for i in 0..length {
            let ref_xor = minimizer[i] ^ key[i];
            let mut mask = 0u8;
            for a in 0..DNA_ALPHABET_SIZE {
                let cand = (a as u8) ^ key[i];
                if cand > ref_xor {
                    mask |= dna_mask(a as u8);
                }
            }
            self.alphabet_i[i] = mask;
        }
        self.alphabet_i[length] = ALL_DNA_MASK;

        self.prefix_letters_vectors[1] = [0; DNA_ALPHABET_SIZE];
        self.prefix_letters_vectors[1][minimizer[0] as usize] = 2;

        for j in 1..length {
            let mut flag = false;
            let mut symb = CmpTag::Eq;

            for i in j..length {
                if !flag {
                    let t = i - j;
                    let left = minimizer[i] ^ key[t];
                    let right = minimizer[t] ^ key[t];
                    if left < right {
                        symb = CmpTag::Lt;
                        flag = true;
                        self.antemer_max_prefix_size = self.antemer_max_prefix_size.min(i + 1);
                        self.postmer_max_size = self.postmer_max_size.min(j + 1);
                    } else if left > right {
                        symb = CmpTag::Gt;
                        flag = true;
                    } else {
                        symb = CmpTag::Eq;
                        let letter = minimizer[t + 1] as usize;
                        self.prefix_letters_vectors[i + 1][letter] =
                            self.prefix_letters_vectors[i + 1][letter].min(j + 1);
                    }
                }
                self.autocorrelation_matrix[i * length + j] = symb;
            }

            for a in 0..DNA_ALPHABET_SIZE {
                if self.prefix_letters_vectors[j + 1][a] == inf {
                    self.prefix_letters_vectors[j + 1][a] =
                        if (a as u8) == minimizer[0] { j + 2 } else { 0 };
                }
            }
        }

        self.suffix_key_convolution
            .resize(self.antemer_max_prefix_size, false);
        self.suffix_key_convolution[..self.antemer_max_prefix_size].fill(false);
        for i in 0..self.antemer_max_prefix_size {
            self.suffix_key_convolution[i] = suffix_key_convolution_gt(minimizer, key, i);
        }

        for i in 1..length {
            let mut correct_alphabet =
                self.alphabet_i[i] & (self.alphabet_i[0] | dna_mask(minimizer[0]));

            for j in 1..i {
                if self.autocorrelation_matrix[(i - 1) * length + j] == CmpTag::Eq {
                    correct_alphabet &= self.alphabet_i[i - j] | dna_mask(minimizer[i - j]);
                }
            }

            for a in 0..DNA_ALPHABET_SIZE {
                let bit = dna_mask(a as u8);
                if correct_alphabet & bit == 0 {
                    continue;
                }
                if self.prefix_letters_vectors[i][a] == 0 {
                    self.antemer_alphabet_zero[i] |= bit;
                } else {
                    self.antemer_alphabet_not_zero[i] |= bit;
                }
            }
        }

        for i in 1..=length {
            for a in 0..DNA_ALPHABET_SIZE {
                let bit = dna_mask(a as u8);
                if self.prefix_letters_vectors[i][a] == 0 {
                    self.postmer_small_values_alphabet_zero[i] |= bit;
                } else {
                    self.postmer_small_values_alphabet_not_zero[i] |= bit;
                }
            }

            if i < length {
                let minimizer_bit = dna_mask(minimizer[i]);
                self.postmer_small_values_alphabet_zero[i] &= !minimizer_bit;
                self.postmer_small_values_alphabet_not_zero[i] &= !minimizer_bit;
            }

            self.postmer_high_values_alphabet[i]
                .push((i, self.alphabet_i[0] | dna_mask(minimizer[0])));
            for k in 1..i {
                if self.autocorrelation_matrix[(i - 1) * length + k] == CmpTag::Eq {
                    self.postmer_high_values_alphabet[i]
                        .push((k, self.alphabet_i[i - k] | dna_mask(minimizer[i - k])));
                }
            }
        }
    }

    pub fn kmer_count(&self, k: usize) -> Result<u128, String> {
        if k < self.length {
            return Err("k must be larger or equal to the length of the minimizer".to_string());
        }

        Ok(self.kmer_count_unchecked(k))
    }

    #[inline]
    pub fn kmer_count_unchecked(&self, k: usize) -> u128 {
        let mut scratch = KmerCountScratch::default();
        self.kmer_count_unchecked_with_scratch(k, &mut scratch)
    }

    #[inline]
    pub fn kmer_count_unchecked_with_scratch(
        &self,
        k: usize,
        scratch: &mut KmerCountScratch,
    ) -> u128 {
        debug_assert!(k >= self.length);

        let alpha = k - self.length;
        let beta_limit = self.postmer_max_size.saturating_sub(2);
        let beta_max = beta_limit.min(alpha);

        self.antemer_with_scratch(alpha, scratch);
        self.postmer_with_scratch(beta_max + self.length, scratch);

        let antemer_array = &scratch.antemer_array[..(alpha + 1)];
        let postmer_array = &scratch.postmer_values[..(beta_max + self.length + 1)];

        let mut total = 0u128;
        for beta in 0..=beta_max {
            total += antemer_array[alpha - beta] * postmer_array[beta + self.length];
        }
        total
    }

    pub fn kmer_counts_up_to(&self, k: usize) -> Result<Vec<u128>, String> {
        if k < self.length {
            return Err("k must be larger or equal to the length of the minimizer".to_string());
        }

        let beta_limit = self.postmer_max_size.saturating_sub(2);
        let antemer_array = self.antemer(k - self.length);
        let postmer_array = self.postmer(beta_limit.min(k - self.length) + self.length);

        let mut values = Vec::with_capacity(k - self.length + 1);
        for k_int in self.length..=k {
            let beta_max_int = beta_limit.min(k_int - self.length);
            let mut total = 0u128;
            for beta in 0..=beta_max_int {
                total +=
                    antemer_array[k_int - self.length - beta] * postmer_array[beta + self.length];
            }
            values.push(total);
        }
        Ok(values)
    }

    fn antemer(&self, alpha: usize) -> Vec<u128> {
        let mut scratch = KmerCountScratch::default();
        self.antemer_with_scratch(alpha, &mut scratch);
        scratch.antemer_array[..(alpha + 1)].to_vec()
    }

    fn antemer_with_scratch(&self, alpha: usize, scratch: &mut KmerCountScratch) {
        let cols = alpha + 1;
        let rows = self.antemer_max_prefix_size;
        prepare_antemer_scratch(scratch, rows, cols);
        let array_prefix = &mut scratch.antemer_prefix[..(rows * cols)];
        let array = &mut scratch.antemer_array[..cols];

        for j in 0..=alpha {
            let mut column_sum = 0u128;
            for i in 0..rows {
                let value = if i > j {
                    0
                } else if j == 0 {
                    1
                } else if i == 0 {
                    dna_mask_len(self.alphabet_i[0]) as u128 * array[j - 1]
                } else if i == j {
                    let mut prod = 1u128;
                    for l in 0..i {
                        let cond = self.autocorrelation_matrix[(i - 1) * self.length + l]
                            == CmpTag::Gt
                            || (self.autocorrelation_matrix[(i - 1) * self.length + l]
                                == CmpTag::Eq
                                && self.suffix_key_convolution[i - 1 - l]);
                        if !cond {
                            prod = 0;
                            break;
                        }
                    }
                    prod
                } else {
                    let mut val =
                        dna_mask_len(self.antemer_alphabet_zero[i]) as u128 * array[j - (i + 1)];

                    for a in 0..DNA_ALPHABET_SIZE {
                        let bit = dna_mask(a as u8);
                        if self.antemer_alphabet_not_zero[i] & bit == 0 {
                            continue;
                        }
                        let t = self.prefix_letters_vectors[i][a];
                        debug_assert!(t > 0 && t <= i + 1);
                        let start = i + 2 - t;
                        let source_col = j - t + 1;
                        for new_prefix in start..rows {
                            val += array_prefix[new_prefix * cols + source_col];
                        }
                    }

                    val
                };
                array_prefix[i * cols + j] = value;
                column_sum += value;
            }
            array[j] = column_sum;
        }
    }

    fn postmer(&self, beta: usize) -> Vec<u128> {
        let mut scratch = KmerCountScratch::default();
        self.postmer_with_scratch(beta, &mut scratch);
        scratch.postmer_values[..(beta + 1)].to_vec()
    }

    fn postmer_with_scratch(&self, beta: usize, scratch: &mut KmerCountScratch) {
        let cols = beta + 1;
        let rows = self.length + 1;
        prepare_postmer_scratch(scratch, rows, cols);
        let array_prefix = &mut scratch.postmer_prefix[..(rows * cols)];
        let array = &mut scratch.postmer_array[..cols];
        let postmer_values = &mut scratch.postmer_values[..cols];

        for j in 0..=beta {
            let mut column_sum = 0u128;
            for i in 0..=self.length {
                let value = if i > j {
                    0
                } else if j == 0 {
                    1
                } else if j < self.length {
                    if i == 0 {
                        (DNA_ALPHABET_SIZE as u128 - 1) * array[j - 1]
                    } else if i == j {
                        1
                    } else {
                        let mut val = dna_mask_len(self.postmer_small_values_alphabet_zero[i])
                            as u128
                            * array[j - (i + 1)];
                        for a in 0..DNA_ALPHABET_SIZE {
                            let bit = dna_mask(a as u8);
                            if self.postmer_small_values_alphabet_not_zero[i] & bit == 0 {
                                continue;
                            }
                            let t = self.prefix_letters_vectors[i][a];
                            debug_assert!(t > 0 && t <= i + 1);
                            let start = i + 2 - t;
                            let end = j + 2 - t;
                            let source_col = j - t + 1;
                            for new_prefix in start..end {
                                val += array_prefix[new_prefix * cols + source_col];
                            }
                        }
                        val
                    }
                } else if j == self.length {
                    if i == self.length {
                        1
                    } else {
                        dna_mask_len(self.alphabet_i[i]) as u128 * pow4(self.length - (i + 1))
                    }
                } else if i == 0 {
                    dna_mask_len(self.alphabet_i[0]) as u128 * array[j - 1]
                } else {
                    let mut correct_alphabet = self.alphabet_i[i];
                    for &(k, mask) in &self.postmer_high_values_alphabet[i] {
                        if j >= self.length + k {
                            correct_alphabet &= mask;
                        }
                    }

                    let mut correct_alphabet_zero = 0u8;
                    let mut correct_alphabet_not_zero = 0u8;

                    for a in 0..DNA_ALPHABET_SIZE {
                        let bit = dna_mask(a as u8);
                        if correct_alphabet & bit == 0 {
                            continue;
                        }
                        let t = self.prefix_letters_vectors[i][a];
                        if t != 0 && j >= self.length - 1 + t {
                            correct_alphabet_not_zero |= bit;
                        } else {
                            correct_alphabet_zero |= bit;
                        }
                    }

                    let mut val = dna_mask_len(correct_alphabet_zero) as u128 * array[j - (i + 1)];

                    for a in 0..DNA_ALPHABET_SIZE {
                        let bit = dna_mask(a as u8);
                        if correct_alphabet_not_zero & bit == 0 {
                            continue;
                        }
                        let t = self.prefix_letters_vectors[i][a];
                        debug_assert!(t > 0);
                        let start = i + 2 - t;
                        let source_col = j - t + 1;
                        for new_prefix in start..=self.length {
                            val += array_prefix[new_prefix * cols + source_col];
                        }
                    }

                    val
                };

                array_prefix[i * cols + j] = value;
                column_sum += value;
            }
            array[j] = column_sum;
            postmer_values[j] = array_prefix[self.length * cols + j];
        }
    }

    pub fn minimizer_len(&self) -> usize {
        self.length
    }
}

#[inline]
fn prepare_antemer_scratch(scratch: &mut KmerCountScratch, rows: usize, cols: usize) {
    let prefix_len = rows * cols;
    if scratch.antemer_prefix.len() < prefix_len {
        scratch.antemer_prefix.resize(prefix_len, 0);
    }
    scratch.antemer_prefix[..prefix_len].fill(0);
    if scratch.antemer_array.len() < cols {
        scratch.antemer_array.resize(cols, 0);
    }
    scratch.antemer_array[..cols].fill(0);
}

#[inline]
fn prepare_postmer_scratch(scratch: &mut KmerCountScratch, rows: usize, cols: usize) {
    let prefix_len = rows * cols;
    if scratch.postmer_prefix.len() < prefix_len {
        scratch.postmer_prefix.resize(prefix_len, 0);
    }
    scratch.postmer_prefix[..prefix_len].fill(0);
    if scratch.postmer_array.len() < cols {
        scratch.postmer_array.resize(cols, 0);
    }
    scratch.postmer_array[..cols].fill(0);
    if scratch.postmer_values.len() < cols {
        scratch.postmer_values.resize(cols, 0);
    }
    scratch.postmer_values[..cols].fill(0);
}

#[inline]
fn suffix_key_convolution_gt(minimizer: &[u8], key: &[u8], i: usize) -> bool {
    let m = minimizer.len();
    if i + 1 >= m {
        return false;
    }

    // Equivalent to:
    // xor_word(minimizer[:m-i-1] + ' ', key[i+1:]) > xor_word(minimizer[i+1:] , key[i+1:])
    // Last character compares ' ' with ' ', so only the first m-i-1 symbols matter.
    let len = m - i - 1;
    for t in 0..len {
        let left = lex_rank(xor_symbol(minimizer[t], key[i + 1 + t]));
        let right = lex_rank(xor_symbol(minimizer[i + 1 + t], key[i + 1 + t]));
        if left > right {
            return true;
        }
        if left < right {
            return false;
        }
    }
    false
}

#[inline]
fn pow4(exp: usize) -> u128 {
    if exp == 0 {
        return 1;
    }
    let shift = exp.checked_mul(2).expect("overflow while computing 4^exp");
    1u128
        .checked_shl(shift as u32)
        .expect("4^exp exceeds u128 capacity")
}
