pub const A: u8 = 0;
pub const C: u8 = 1;
pub const G: u8 = 2;
pub const T: u8 = 3;
pub const SPACE: u8 = 4;

pub const DNA_ALPHABET_SIZE: usize = 4;
pub const ALL_DNA_MASK: u8 = 0b1111;

#[inline]
pub const fn dna_mask(symbol: u8) -> u8 {
    1u8 << symbol
}

#[inline]
pub const fn dna_mask_len(mask: u8) -> u32 {
    mask.count_ones()
}

#[inline]
pub fn xor_symbol(left: u8, right: u8) -> u8 {
    if left == SPACE || right == SPACE {
        SPACE
    } else {
        left ^ right
    }
}

#[inline]
pub const fn lex_rank(symbol: u8) -> u8 {
    // Keeps the same order as Python string comparison: ' ' < 'A' < 'C' < 'G' < 'T'.
    if symbol == SPACE { 0 } else { symbol + 1 }
}

#[inline]
pub fn symbol_from_char(c: char) -> Option<u8> {
    match c {
        'A' => Some(A),
        'C' => Some(C),
        'G' => Some(G),
        'T' => Some(T),
        _ => None,
    }
}

#[inline]
pub const fn char_from_symbol(symbol: u8) -> char {
    match symbol {
        A => 'A',
        C => 'C',
        G => 'G',
        T => 'T',
        _ => ' ',
    }
}

pub fn parse_dna_word(word: &str) -> Result<Vec<u8>, String> {
    if word.is_empty() {
        return Err("DNA word cannot be empty".to_string());
    }

    let mut out = Vec::with_capacity(word.len());
    for c in word.chars() {
        let symbol = symbol_from_char(c).ok_or_else(|| {
            format!("Invalid DNA symbol '{c}'. Allowed symbols are A, C, G, T (uppercase).")
        })?;
        out.push(symbol);
    }
    Ok(out)
}

pub fn format_dna_word(symbols: &[u8]) -> String {
    symbols.iter().map(|&s| char_from_symbol(s)).collect()
}

pub fn decode_index_to_kmer(mut index: u64, m: usize) -> Vec<u8> {
    let mut out = vec![A; m];
    for i in (0..m).rev() {
        out[i] = (index & 0b11) as u8;
        index >>= 2;
    }
    out
}

#[inline]
pub fn decode_index_to_kmer_inplace(mut index: u64, out: &mut [u8]) {
    for i in (0..out.len()).rev() {
        out[i] = (index & 0b11) as u8;
        index >>= 2;
    }
}

pub fn total_kmers(m: usize) -> Result<u64, String> {
    if m > 31 {
        return Err(
            "m is too large for full enumeration in this build (m must be <= 31 for 64-bit indexing)"
                .to_string(),
        );
    }
    Ok(1u64 << (2 * m))
}
