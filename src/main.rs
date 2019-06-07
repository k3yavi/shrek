
fn main() {
    let matches = App::new("debruijn graph")
        .version("1.0")
        .author("Avi S. <asrivastava@cs.stonybrook.edu>")
        .about("De-bruijn graph creation")
        .arg(
            Arg::with_name("fasta")
                .short("f")
                .long("fasta")
                .value_name("FILE")
                .help(
                    "Txome/Genome Input Fasta file, (Needed only with -m i.e. while making index)",
                ),
        ).arg(
            Arg::with_name("index")
                .short("i")
                .long("index")
                .value_name("FILE")
                .help("Index of the reference")
                .required(true),
        ).get_matches();

    // initializing logger
    pretty_env_logger::init_timed();

    // obtain reader or fail with error (via the unwrap method)
    let index_file = matches.values_of("index").unwrap().next().unwrap();

    // Gets a value for config if supplied by user
    let fasta_file = matches.value_of("fasta").unwrap();
    info!("Path for reference FASTA: {}", fasta_file);

    // if index not found then create a new one
    let reader = fasta::Reader::from_file(fasta_file).unwrap();
    let (seqs, tgmap, gene_order) = read_fasta(reader, tgmap_file);

    //Set up the filter_kmer call based on the number of sequences.
    filter_kmers_callback(&seqs, index_file, &uhs, tgmap, gene_order);

    info!("Finished Indexing !");
}
