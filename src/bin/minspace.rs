use std;
use seq_io;
use seq_io::prelude::*;
use nthash;

fn encode_rle(inp_seq: &str) -> (String, Vec<usize>) {
    let mut prev_char = '#';
    let mut hpc_seq = String::new();
    let mut pos_vec = Vec::<usize>::new();
    let mut prev_i = 0;
    for (i, c) in inp_seq.chars().enumerate() {
        if c == prev_char && "ACTGactgNn".contains(c) {continue;}
        if prev_char != '#' {
            hpc_seq.push(prev_char);
            pos_vec.push(prev_i);
            prev_i = i;
        }
        prev_char = c;
    }
    hpc_seq.push(prev_char);
    pos_vec.push(prev_i);
    (hpc_seq, pos_vec)
}

pub fn lmers_hashes(inp_seq_raw: String, l: usize, delta: f64) -> Vec<(String, u64)> {
    let mut result = Vec::<(String, u64)>::new();
    let hash_bound = (delta * (u64::max_value() as f64)) as u64;
    let tup = encode_rle(&inp_seq_raw);
    let inp_seq = tup.0;
    let minimizers_with_pos = nthash::NtHashIterator::new(inp_seq.as_bytes(), l)
        .unwrap().enumerate().filter(|(_, x)| *x < hash_bound);
    for (i, hash) in minimizers_with_pos {
        let lmer = String::from(&inp_seq[i..i+l]);
        result.push((lmer, hash));
    }
    result
}

pub fn kminmers(lmers_hashes: Vec<(String, u64)>, k: usize) -> Vec<Vec<u64>> {
    if lmers_hashes.len() < k {
        vec![lmers_hashes.iter().map(|p| p.1).collect::<Vec<_>>()]
    } else {
        lmers_hashes
            .windows(k)
            .map(|w| w.iter().map(|p| p.1).collect::<Vec<_>>())
            .collect()
    }
}


fn main() {
    let mut reader = seq_io::fasta::Reader::from_path("/home/reepoi/GitHub/rust-mdbg/example/reads-0.00.fa").unwrap();
    while let Some(record) = reader.next() {
        let record = record.expect("Error reading record");
        let seq = std::str::from_utf8(record.seq()).unwrap().to_string();
        let result = lmers_hashes(seq, 10, 0.0001);
        println!("lmers and hashes: {:?}", result);
        println!("kminmers: {:?}", kminmers(result, 4));
        break;
    }
}
