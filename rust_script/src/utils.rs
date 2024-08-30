use std::collections::HashMap;
use num::One;
use num::Zero;
use rand::Rng;
use num::pow;
use num::BigUint;
use num::range_inclusive;

pub fn binary_vec_to_kmer(bin_vec:Vec<u64>,k:usize)-> String {

    assert!(bin_vec.len()==2*k);

    let mut s = String::new();

    for i in 0..k {

        let bits = (bin_vec[2*i],bin_vec[2*i+1]);

        match bits {
            (0,0) => s.push('A'),
            (0,1) => s.push('C'),
            (1,0)=> s.push('G'),
            (1,1)=> s.push('T'),
            _=> panic!()
        }

    }

    assert!(s.len()==k);

    return s;


}


pub fn generate_random_k_mers(k:usize, n:usize) -> Vec<String> {

    // let nmax: BigUint = pow(BigUint::from(2u32), 2*k)-BigUint::one();

    let mut l: Vec<String> = vec![];

    let mut rng = rand::thread_rng();

    for _ in 0..n {
        let bin_vec: Vec<u64> = (0..2*k).map(|_| rng.gen_range(0, 2)).collect();
        l.push(binary_vec_to_kmer(bin_vec, k));
    }
    return l;
}

pub fn biguint_to_binary_vec(n: BigUint,k:usize)-> Vec<u64> {

    let mut bin_vec = vec![];

    let x = &mut n.clone();

    while x!= &BigUint::zero() {

        let y = x.clone();

        let rem = y % BigUint::from(2u64);
        let bit = rem.to_u64_digits();

        if bit == vec![] {
            bin_vec.insert(0,0);
        }
        else {
            bin_vec.insert(0,1);
        }

        *x /= BigUint::from(2u64);
    }


    while bin_vec.len()<2*k {
        bin_vec.insert(0,0);
    }

    return bin_vec

}

pub fn generate_all_k_mers(k:usize) -> Vec<String> {

    let mut l: Vec<String> = vec![];

    let max = pow(BigUint::from(4u8),k) - BigUint::one();

    for i in range_inclusive(BigUint::from(0u64),max) {

        l.push(binary_vec_to_kmer(biguint_to_binary_vec(i,k), k));

    }

    return l;
}

pub fn find_minimizer(s:&String, m: usize) -> usize {

    assert!(s.len()>=m, "String too short");

    let mut index_min = 0;

    for i in 1..(s.len()-m+1) {
        for j in 0..m {

            let a = s.chars().nth(i+j);
            let b = s.chars().nth(index_min+j);

            if a < b {
                index_min = i;
                break;
            }
            else if b <  a {
                break;
            }
        }
    }

    return index_min;
}

pub fn number_of_greater_letters(alphabet:&Vec<char>) -> HashMap<char, usize> {

    let mut dic: HashMap<char, usize> = HashMap::new();

    for a in alphabet.iter() {

        let mut count: usize = 0;

        for b in alphabet.iter(){
            if b>a {
                count +=1;
            }
        }

        dic.insert(a.to_ascii_uppercase(), count);
    }

    return dic;
}


pub fn number_of_greater_words(s:&String,alphabet:&Vec<char>) -> BigUint {
    // greater or equal

    let greater_letters_dic = number_of_greater_letters(alphabet);

    let m = s.len();
    let mut sum = BigUint::one(); // pour le préfixe=m

    for (i, c) in s.chars().enumerate() {
        let prod = greater_letters_dic.get(&c).unwrap() * pow(BigUint::from(alphabet.len()),m-i-1);
        sum+= prod;
    }

    return sum;
}