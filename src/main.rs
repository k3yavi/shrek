extern crate bio;
extern crate clap;
extern crate debruijn;
extern crate boomphf;
extern crate pretty_env_logger;

#[macro_use]
extern crate log;

#[macro_use]
extern crate lazy_static;

use std::io;
use std::fs::{OpenOptions};
use std::io::Write;
use std::sync::Arc;
use bio::io::{fasta};
use clap::{App, Arg, ArgMatches, SubCommand};

use debruijn::*;
use debruijn::filter::*;
use debruijn::graph::*;
use debruijn::compression::*;
use debruijn::dna_string::{DnaString, DnaStringSlice};

use boomphf::hashmap::{BoomHashMap2, NoKeyBoomHashMap};
use rayon::prelude::*;

mod utils;

const MIN_KMERS: usize = 1;
const MAX_WORKER: usize = 10;
pub const MEM_SIZE: usize = 1;
pub const STRANDED: bool = false;
const U32_MAX: u32 = u32::max_value();
pub const REPORT_ALL_KMER: bool = false;
const MIN_SHARD_SEQUENCES: usize = 0;

type KmerType = kmer::Kmer31;
type PmerType = debruijn::kmer::Kmer6;
lazy_static! {
    static ref PERM: Vec<usize> = {
        let maxp = 1 << (2 * PmerType::k());
        let mut permutation = Vec::with_capacity(maxp);
        for i in 0..maxp {
            permutation.push(i);
        }
        permutation
    };
}

fn group_by_slices<T, K: PartialEq, F: Fn(&T) -> K>(
    data: &[T],
    f: F,
    min_size: usize,
) -> Vec<&[T]> {
    let mut slice_start = 0;
    let mut result = Vec::new();
    for i in 1..data.len() {
        if !(f(&data[i - 1]) == f(&data[i])) && (i - slice_start) > min_size {
            result.push(&data[slice_start..i]);
            slice_start = i;
        }
    }

    if slice_start > 0 || (data.len() - slice_start) > min_size {
        result.push(&data[slice_start..]);
    }
    result
}

fn assemble_shard<K: Kmer>(
    shard_data: &[(u16, u32, DnaStringSlice, Exts)],
    summarizer: &Arc<CountFilterEqClass<u32>>,
) -> BaseGraph<K, (EqClassIdType, u8)> {
    let filter_input: Vec<_> = shard_data
        .into_iter()
        .cloned()
        .map(|(_, seqid, string, exts)| (string, exts, seqid))
        .collect();

    let (phf, _): (BoomHashMap2<K, Exts, (EqClassIdType, u8)>, _) = filter_kmers(
        &filter_input,
        summarizer,
        STRANDED,
        REPORT_ALL_KMER,
        MEM_SIZE,
    );

    //println!("printing filters");
    //println!("{:?}", phf);
    compress_kmers_with_hash(STRANDED, ScmapCompress::new(), &phf)
}

fn merge_shard_dbgs<K: Kmer + Sync + Send>(
    uncompressed_dbgs: Vec<BaseGraph<K, (EqClassIdType, u8)>>,
) -> DebruijnGraph<K, (EqClassIdType, u8)> {
    let combined_graph = BaseGraph::combine(uncompressed_dbgs.into_iter()).finish();
    compress_graph(STRANDED, ScmapCompress::new(), combined_graph, None)
}

fn partition_contigs<'a, K: Kmer>(
    contig: &'a DnaString,
    contig_id: u32,
) -> Vec<(u16, u32, DnaStringSlice<'a>, Exts)> {
    // One FASTA entry possibly broken into multiple contigs
    // based on the location of `N` int he sequence.

    let mut bucket_slices = Vec::new();

    if contig.len() >= K::k() {
        let msps = debruijn::msp::simple_scan::<_, PmerType>(K::k(), contig, &PERM, !STRANDED);
        for msp in msps {
            let bucket_id = msp.bucket();
            let slice = contig.slice(msp.start(), msp.end());
            let exts = Exts::from_dna_string(contig, msp.start(), msp.len());
            bucket_slices.push((bucket_id, contig_id, slice, exts));
        }
    }

    bucket_slices
}

fn generate(sub_m: &ArgMatches) -> Result<(), io::Error> {
    // obtain reader or fail with error (via the unwrap method)
    let gfa_file = sub_m.values_of("gfa").unwrap().next().unwrap();
    info!("GFA will be created at: {}", gfa_file);

    // Gets a value for config if supplied by user
    let fasta_file = sub_m.value_of("fasta").unwrap();
    info!("Path for reference FASTA: {}", fasta_file);

    let reader = fasta::Reader::from_file(fasta_file).unwrap();

    let mut seqs = Vec::new();
    let mut seq_names = Vec::new();
    {
        // scope for reading data
        let mut transcript_counter = 0;

        info!("Starting reading the Fasta file\n");
        for result in reader.records() {
            // obtain record or fail with error
            let record = result.unwrap();

            // Sequence
            let dna_string = DnaString::from_acgt_bytes_hashn(record.seq(), record.id().as_bytes ());
            seqs.push(dna_string);
            seq_names.push( record.id().to_string() );

            transcript_counter += 1;
            if transcript_counter % 100 == 0 {
                print!("\r Done reading {} sequences", transcript_counter);
                io::stdout().flush().expect("Could not flush stdout");
            }
        }

        println!();
        info!(
            "Done reading the Fasta file; Found {} sequences",
            transcript_counter
        );
    }

    if seqs.len() >= std::u32::MAX as usize {
        panic!("Too many ({}) sequences to handle.", seqs.len());
    }

    info!("Sharding sequences...");
    let mut buckets: Vec<_> = seqs
        .par_iter()
        .enumerate()
        .flat_map(|(id, seq)| partition_contigs::<KmerType>(seq, id as u32))
        .collect();

    buckets.par_sort_unstable_by_key(|x| x.0);
    info!("Got {} sequence chunks", buckets.len());

    let summarizer = Arc::new(debruijn::filter::CountFilterEqClass::new(MIN_KMERS));
    let sequence_shards = group_by_slices(&buckets, |x| x.0, MIN_SHARD_SEQUENCES);
    let mut shard_dbgs = Vec::with_capacity(sequence_shards.len());

    info!("Assembling {} shards...", sequence_shards.len());
    sequence_shards
        .into_par_iter()
        .map_with(summarizer.clone(), |s, strings| {
            assemble_shard::<KmerType>(strings, s)
        }).collect_into_vec(&mut shard_dbgs);

    println!();
    //println!("{:?}", shard_dbgs);
    info!("Done separate de Bruijn graph construction");
    info!("Starting merging disjoint graphs");
    let dbg = merge_shard_dbgs(shard_dbgs);
    let eq_classes = summarizer.get_eq_classes();

    info!("Graph merge complete");
    info!("Writing GFA !");
    write_gfa(eq_classes, dbg, gfa_file, seqs, seq_names)
        .expect("Can't write gfa");

    Ok(())
}

pub fn write_gfa( _eq_classes: Vec<Vec<u32>>,
                  dbg: DebruijnGraph<KmerType, (u32, u8)>,
                  gfa_file: &str,
                  seqs: Vec<DnaString>,
                  seq_names: Vec<String>)
                  -> Result<(), io::Error> {
    // writing S and L flags into the gfa
    dbg.to_gfa(gfa_file)?;
    info!("GFA: S and L tag written !");

    let dbg_index;
    {
        let mut total_kmers = 0;
        let kmer_length = KmerType::k();
        for node in dbg.iter_nodes() {
            total_kmers += node.len() - kmer_length + 1;
        }

        info!("Total {:?} kmers to process in dbg", total_kmers);
        info!("Making mphf of kmers");
        let mphf = boomphf::Mphf::from_chunked_iterator_parallel(1.7, &dbg, None, total_kmers, MAX_WORKER);

        info!("Assigning offsets to kmers");
        let mut node_and_offsets = Vec::with_capacity(total_kmers);
        node_and_offsets.resize(total_kmers, (U32_MAX as u32, U32_MAX as u32));

        for node in &dbg {
            let node_id = node.node_id;

            for (offset, kmer) in node.into_iter().enumerate() {
                let index = match mphf.try_hash(&kmer) {
                    None => panic!("can't find kmer"),
                    Some(index) => index,
                };

                node_and_offsets[index as usize] = (node_id as u32, offset as u32);
            }
        }

        info!("Done creating index");
        dbg_index = NoKeyBoomHashMap::new_with_mphf(mphf, node_and_offsets);
    }

    let mut path: String = "".to_string();
    let mut wtr = OpenOptions::new().append(true).open(gfa_file)?;
    for (seq_index, seq) in seqs.into_iter().enumerate() {
        let mut coverage: usize = 0;
        let mut num_nodes: usize = 0;
        let mut kmer_it = seq.iter_kmers::<KmerType>();

        path.clear();
        let mut old_kmer = None;
        while let Some(kmer) = kmer_it.next() {
            let mut is_kmer_rc = false;
            let kmer_rc = kmer.rc();

            let fwd_idx = match dbg_index.get(&kmer) {
                Some((nid, offset)) => {
                    let node = dbg.get_node(*nid as usize);
                    let seq = node
                        .sequence()
                        .get_kmer(*offset as usize);
                    if kmer == seq { Some(node) } else { None }
                },
                None => None,
            };

            let rev_idx = match dbg_index.get(&kmer_rc) {
                Some((nid, offset)) => {
                    let node = dbg.get_node(*nid as usize);
                    let seq = node
                        .sequence()
                        .rc()
                        .get_kmer(*offset as usize);
                    if kmer_rc == seq { is_kmer_rc = true; Some(node) } else { None }
                },
                None => None,
            };

            if fwd_idx.is_some() && rev_idx.is_some() {
                println!("fwd: {:?} {:?}", fwd_idx, kmer);
                println!("rev: {:?} {:?}", rev_idx, kmer_rc);
                panic!("found both fwd and rev kmer");
            }

            if fwd_idx.is_none() && rev_idx.is_none() {
                panic!("Neither fwd nor reverse kmer found");
            }

            let node = if fwd_idx.is_some() { fwd_idx.unwrap() } else {rev_idx.unwrap()};
            let node_len = node.len();
            let node_sign = match is_kmer_rc {
                true => "-",
                false => "+",
            };

            num_nodes += 1;
            match num_nodes {
                1 => coverage += node_len,
                _ => coverage += node_len - KmerType::k() + 1,
            }

            if path.is_empty() {
                path = format!("{}{}", node.node_id, node_sign);
            } else {
                path = format!("{},{}{}", path, node.node_id, node_sign);
            }

            let mut skip_kmers = node_len - KmerType::k() ;
            while skip_kmers > 0 {
                kmer_it.next();
                skip_kmers -= 1;
            }

            if old_kmer.is_some() && old_kmer == Some(kmer) { break; }
            old_kmer = Some(kmer);
            assert!(coverage <= seq.len(), "{} {} {} {} {:?}, {:?}, {:?} {:?}",
                    num_nodes, coverage, seq.len(),
                    node_len - KmerType::k() + 1, kmer,
                    node, node.sequence(), seq);
        } // end-while

        let seq_name = seq_names.get(seq_index).unwrap();
        //assert!(coverage == seq.len(),
        //        "didn't cover full transcript {}: len: {} covered: {}",
        //        seq_name, seq.len(), coverage);

        writeln!(wtr, "P\t{}\t{}\t*",
                 seq_name,
                 path
        )?;
    }// end- seq for
    info!("Done writing Path");

    Ok(())
}

fn main() -> io::Result<()> {
    let matches = App::new("Shrek")
        .version("1.0")
        .author("Avi S. <asrivastava@cs.stonybrook.edu>")
        .about("oh this is another one of those onion things !")
        .subcommand(
            SubCommand::with_name("generate")
                .arg(
                    Arg::with_name("fasta")
                        .short("f")
                        .long("fasta")
                        .value_name("FILE")
                        .help("Txome/Genome Input fasta")
                        .required (true),
                ).arg(
                    Arg::with_name("gfa")
                        .short("g")
                        .long("gfa")
                        .value_name("FILE")
                        .help("path to output gfa file")
                        .required(true),
                ),
        )
        .subcommand(
            SubCommand::with_name("compare")
                .arg(
                    Arg::with_name("gfa1")
                        .short("g1")
                        .long("gfa1")
                        .value_name("FILE")
                        .help("GFA1 to compare")
                        .required (true),
                ).arg(
                    Arg::with_name("gfa2")
                        .short("g2")
                        .long("gfa2")
                        .value_name("FILE")
                        .help("GFA2")
                        .required(true),
                ),
        ).get_matches();

    // initializing logger
    pretty_env_logger::init_timed();

    match matches.subcommand_matches("generate") {
        Some(sub_m) => {
            let ret = generate(&sub_m);
            return ret;
        }
        None => (),
    };

    match matches.subcommand_matches("compare") {
        Some(sub_m) => {
            let ret = utils::compare(&sub_m);
            return ret;
        }
        None => (),
    };

    Ok(())
}
