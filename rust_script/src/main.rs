// use std::io;

use substring::Substring;
mod utils;

fn main() {

    let n = 5;
    let k = 5;
    let m = 3;

    let kmers_list = utils::generate_random_k_mers(k,n);

    for kmer in kmers_list.iter() {

        let min_index = utils::find_minimizer(kmer, m);

        let minimizer = kmer.substring(min_index,min_index+m);

        println!("{kmer} : {minimizer}");
    }

    // println!("=====");

    // let kmers_list = utils::generate_all_k_mers(k);

    // for kmer in kmers_list.iter() {

    //     println!("{kmer}");
    // }

    
}

