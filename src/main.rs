use std::ffi::CStr;

use hmmer_rs::*;

use log::*;
use env_logger;

fn main() {
    env_logger::init();

    let hmm = Hmm::read_one_from_path(std::path::Path::new("tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm")).unwrap();

    println!("HMM name: {}", hmm.name());

    let mut hmmsearch = HmmerPipeline::new(&hmm);

    let mut query_seq = EaselSequence::new(hmm.c_alphabet());
    let seq: &[u8] = "MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ".as_bytes();

    query_seq.replace_sequence(&seq);
    debug!("Query seq replaced;");
    error!("It appears that setting the sequence like does not work, so no results are returned. Need to bugfix.");
    hmmsearch.query(&query_seq);

    let hmmsearch_result = hmmsearch.get_results();

    debug!("HMMsearch result: {:?}", hmmsearch_result);

    // th->hit[h]->score,
    println!("First domain nreported: {}", unsafe {(*hmmsearch_result.c_th).nreported});
    let first_hit = unsafe {*(*(*hmmsearch_result.c_th).hit.offset(0))};
    println!("Name of first hit {}", unsafe { CStr::from_ptr(first_hit.name).to_string_lossy() });
    println!("Score of first hit overall {}", first_hit.score);

    let first_domain = unsafe {*first_hit.dcl.offset(0)};
    println!("First domain score: {}", first_domain.bitscore);
    // exp(th->hit[h]->lnP) * pli->Z;
    let evalue = first_domain.lnP.exp() * unsafe {(*hmmsearch_result.c_pli).Z};
    println!("First domain evalue: {:?}", evalue);
}

