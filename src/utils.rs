use std::io;
use std::fs::File;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};

use clap::{ArgMatches};

pub fn get_data(gfa_file: &str)
                -> (HashMap<usize, usize>, HashMap<String, String>){
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
                let rseq: usize = toks[2].len();
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

pub fn get_lens(u1: &HashMap<usize, usize>, path_seq: &String)
                -> usize {
    let mut total_len = 0;
    for path in path_seq.trim().split(",") {
        let rid: usize = path
            .trim_matches(|c| c == '-' || c == '+')
            .parse::<usize>()
            .unwrap();

        total_len += u1.get(&rid).unwrap();
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

    for (id1, path1) in p1 {
        let len1 = get_lens(&u1, &path1);

        let path2 = p2.get(&id1).unwrap();
        let len2 = get_lens(&u2, &path2);

        assert!(len1 == len2);
    }

    info!("All Done!!");
    Ok(())
}
