use num::One;
use num::Zero;
use rand::seq::index;
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


// def number_of_greater_letters(alphabet = { 'A', 'T', 'C', 'G'}):
//     dic = {}
//     for a in alphabet:
//         dic[a] = len([s for s in alphabet if s>a])
//     return dic

// def number_of_greater_words(string:str,alphabet = { 'A', 'T', 'C', 'G'}):
//     # greater or equal
//     greater_letters_dic = number_of_greater_letters(alphabet)

//     m  = len(string)
//     sum = 1 # pour le préfixe=m
//     for i in range(0,m):
//         ai = string[i]
//         prod = greater_letters_dic[ai]*len(alphabet)**(m-i-1)
//         sum+= prod
//     return sum