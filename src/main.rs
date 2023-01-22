#![allow(unused_imports)]

use std::ffi::CStr;

use hmmer_rs::*;

use env_logger;
use log::*;

fn main() {
    env_logger::init();
    let hmms = Hmm::read_hmms_from_path(std::path::Path::new(
        "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm",
    )).unwrap();
    let hmm = &hmms[0];

    println!("HMM name: {}", hmm.name());

    let mut hmmsearch = HmmerPipeline::new(&hmm);

    let mut query_seq = EaselSequence::new(hmm.c_alphabet());
    let seq: &[u8] =
        b"MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ";

    query_seq.replace_sequence(&seq).unwrap();
    debug!("Query seq replaced;");

    hmmsearch.query(&query_seq);

    let hmmsearch_result = hmmsearch.get_results();

    println!("Total number of reported hits: {}", hmmsearch_result.nreported());

    for hit in hmmsearch_result.hits() {
        println!("New hit:");
        println!("Hit name: {}", hit.name());
        println!("Hit score: {}", hit.score());
        for domain in hit {
            println!("New domain:");
            println!("Domain score: {}", domain.bitscore());
            println!("Domain evalue: {:?}", domain.evalue());
        }
    }
}
