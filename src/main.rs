use hmmer_rs::*;

use log::*;
use env_logger;

use std::ffi::CStr;

fn main() {
    env_logger::init();
    
    let hmm = Hmm::read_one_from_path(std::path::Path::new("test/data/DNGNGWU00010_mingle_output_good_seqs.hmm")).unwrap();

    println!("HMM name: {}", unsafe { CStr::from_ptr((*hmm.c_hmm).name).to_string_lossy() });

    let hmmsearch = HmmerPipeline{};

    let hmmsearch_result = hmmsearch.run_hmm_on_file(&hmm, std::path::Path::new("test/data/graftm4o5_y58f.head2.faa"));

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

