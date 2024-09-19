// use std::io;

use substring::Substring;
mod utils;
mod lib;

fn main() {

    let alphabet = vec!['A','C','G','T'];
    let n = 1;
    let k = 11;
    let m = 6;

    let kmers_list = utils::generate_random_k_mers(k,n);
    // let kmers_list = utils::generate_all_k_mers(k);

    // let test= utils::number_of_greater_words(&"ACACAC".to_string(),&alphabet);

    // println!("{test}");

    for kmer in kmers_list.iter() {

        let min_index = utils::find_minimizer(kmer, m);

        let minimizer = &kmer.substring(min_index,min_index+m).to_string();

        let g = utils::number_of_greater_words(minimizer, &alphabet);

        println!("{kmer} : {minimizer}  {g} ");

        lib::autocorrelation_matrix(minimizer, &alphabet);
    }
    
}

