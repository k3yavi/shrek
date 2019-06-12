use std::io;
use std::fs::File;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};

use debruijn::Mer;
use clap::{ArgMatches};

use debruijn::dna_string::{DnaString};

pub fn get_data(gfa_file: &str)
                -> (HashMap<usize, String>, HashMap<String, String>){
    let f = BufReader::new(File::open(gfa_file).unwrap());
    let mut unitigs = HashMap::new();
    let mut path = HashMap::new();

    for line in f.lines() {
        let toks: Vec<String> = line.expect("Unable to read line")
            .trim()
            .split("\t")
            .map(|chunk| chunk.to_string())
            .collect();

        let rtype = &toks[0];
        match rtype.as_str() {
            "S" => {
                let rid: usize = toks[1].parse().unwrap();
                let rseq: String = toks[2].clone();
                unitigs.insert(rid, rseq);
            }
            "P" => {
                let rname = toks[1].clone();
                let rpath = toks[2].clone();
                path.insert(rname, rpath);
            }
            _ => ()
        }
    }

    return (unitigs, path)
}

pub fn get_lens(u1: &HashMap<usize, String>, path_seq: &String)
                -> Vec<String> {
    let mut total_len = Vec::new();
    for path in path_seq.trim().split(",") {
        let rid: usize = path
            .trim_matches(|c| c == '-' || c == '+')
            .parse::<usize>()
            .unwrap();

        total_len.push( u1.get(&rid).unwrap().clone() );
    }

    total_len
}

pub fn compare(sub_m: &ArgMatches) -> Result<(), io::Error> {
    // obtain reader or fail with error (via the unwrap method)
    let gfa1_file = sub_m.values_of("gfa1").unwrap()
        .next().unwrap();
    let gfa2_file = sub_m.values_of("gfa2").unwrap()
        .next().unwrap();

    info!("Parsing {}", gfa1_file);
    let (u1, p1) = get_data(gfa1_file);

    info!("Parsing {}", gfa2_file);
    let (u2, p2) = get_data(gfa2_file);

    //assert!(u1.len() == u2.len(), "{:?}, {:?}", u1.len(), u2.len());
    let mut diff = 0;
    let p1_len = p1.len();
    for (id1, path1) in p1 {
        let len1 = get_lens(&u1, &path1);

        let path2 = p2.get(&id1).unwrap();
        let len2 = get_lens(&u2, &path2);

        let mut total_length1: usize = 0;
        let mut total_length2: usize = 0;
        for unitig in len1 { total_length1 += unitig.len(); }
        for unitig in len2 { total_length2 += unitig.len(); }

        if total_length1 != total_length2 {
            diff += 1;
            println!("{:?}, {:?} {:?}", id1, total_length1, total_length2);
        }
        //if len1.len() != len2.len() {
        //    println!("Overall: {:?} \n {:?} \n {:?}", id1, len1.len(), len2.len());
        //    //println!("Overall: {:?} \n {:?} \n {:?}", id1, len1, len2);
        //    diff += 1;
        //    //break;
        //} else {
        //    for i in 0..len1.len() {
        //        let seq1 = DnaString::from_acgt_bytes_hashn(len1[i].as_bytes(), &[i as u8]);
        //        let seq2 = DnaString::from_acgt_bytes_hashn(len2[i].as_bytes(), &[i as u8]);
        //        if seq1 != seq2 &&
        //            seq1.rc() != seq2.rc() &&
        //            seq1 != seq2.rc() &&
        //            seq1.rc() != seq2 {
        //                diff += 1;
        //                println!("{:?} {}", seq1, len1[i].len());
        //                println!("{:?} {}", seq2.rc(), len2[i].len());
        //        }
        //    }
        //}
    }

    info!("All Done!! Diff = {}/{}", diff, p1_len);
    Ok(())
}
