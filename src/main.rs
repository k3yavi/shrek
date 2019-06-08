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
use std::io::Write;
use std::sync::Arc;
use bio::io::{fasta};
use clap::{App, Arg, ArgMatches, SubCommand};

use debruijn::*;
use debruijn::filter::*;
use debruijn::graph::*;
use debruijn::compression::*;
use debruijn::dna_string::{DnaString, DnaStringSlice};

use boomphf::hashmap::BoomHashMap2;
use rayon::prelude::*;

mod utils;

const MIN_KMERS: usize = 1;
pub const MEM_SIZE: usize = 1;
pub const STRANDED: bool = true;
pub const REPORT_ALL_KMER: bool = false;
const MIN_SHARD_SEQUENCES: usize = 2000;

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
    if slice_start > 0 {
        result.push(&data[slice_start..]);
    }
    result
}

fn assemble_shard<K: Kmer>(
    shard_data: &[(u16, u32, DnaStringSlice, Exts)],
    summarizer: &Arc<CountFilterEqClass<u32>>,
) -> BaseGraph<K, EqClassIdType> {
    let filter_input: Vec<_> = shard_data
        .into_iter()
        .cloned()
        .map(|(_, seqid, string, exts)| (string, exts, seqid))
        .collect();

    let (phf, _): (BoomHashMap2<K, Exts, EqClassIdType>, _) = filter_kmers(
        &filter_input,
        summarizer,
        STRANDED,
        REPORT_ALL_KMER,
        MEM_SIZE,
    );

    compress_kmers_with_hash(STRANDED, ScmapCompress::new(), &phf)
}

fn merge_shard_dbgs<K: Kmer + Sync + Send>(
    uncompressed_dbgs: Vec<BaseGraph<K, EqClassIdType>>,
) -> DebruijnGraph<K, EqClassIdType> {
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
        let msps = debruijn::msp::simple_scan::<_, PmerType>(K::k(), contig, &PERM, false);
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

    info!("Done separate de Bruijn graph construction");
    info!("Starting merging disjoint graphs");
    let dbg = merge_shard_dbgs(shard_dbgs);

    info!("Graph merge complete");
    dbg.to_gfa(gfa_file).expect("can't write gfa");

    info!("Finished Indexing !");

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
