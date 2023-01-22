use std::ffi::CStr;

use hmmer_rs::*;
use log::*;

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn hmmsearch_on_file() {
        let hmms = Hmm::read_hmms_from_path(std::path::Path::new(
            "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm",
        )).unwrap();
        let hmm = &hmms[0];
        
        println!("HMM name: {}", unsafe {
            CStr::from_ptr((*hmm.c_hmm).name).to_string_lossy()
        });

        let mut hmmsearch = HmmerPipeline::new(&hmm);

        let hmmsearch_result = hmmsearch.run_hmm_on_file(
            &hmm,
            std::path::Path::new("tests/data/graftm4o5_y58f.head2.faa"),
        );

        debug!("HMMsearch result: {:?}", hmmsearch_result);

        // th->hit[h]->score,
        println!("First domain nreported: {}", unsafe {
            (*hmmsearch_result.c_th).nreported
        });
        let first_hit = unsafe { *(*(*hmmsearch_result.c_th).hit.offset(0)) };
        println!("Name of first hit {}", unsafe {
            CStr::from_ptr(first_hit.name).to_string_lossy()
        });
        println!("Score of first hit overall {}", first_hit.score);

        let first_domain = unsafe { *first_hit.dcl.offset(0) };
        println!("First domain score: {}", first_domain.bitscore);
        // exp(th->hit[h]->lnP) * pli->Z;
        let evalue = first_domain.lnP.exp() * unsafe { (*hmmsearch_result.c_pli).Z };
        println!("First domain evalue: {:?}", evalue);

        assert_eq!(hmmsearch_result.nreported(), 1);
        assert_eq!(evalue, 1.4970530541655288e-48);
    }

    #[test]
    fn test_hmmsearch_by_query() {
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
        error!("It appears that setting the sequence like does not work, so no results are returned. Need to bugfix.");

        hmmsearch.query(&query_seq);

        let hmmsearch_result = hmmsearch.get_results();

        debug!("HMMsearch result: {:?}", hmmsearch_result);

        // th->hit[h]->score,
        println!("nreported: {}", hmmsearch_result.nreported());
        let first_hit = unsafe { *(*(*hmmsearch_result.c_th).hit.offset(0)) };
        println!("Name of first hit {}", unsafe {
            CStr::from_ptr(first_hit.name).to_string_lossy()
        });
        println!("Score of first hit overall {}", first_hit.score);

        let first_domain = unsafe { *first_hit.dcl.offset(0) };
        println!("First domain score: {}", first_domain.bitscore);
        // exp(th->hit[h]->lnP) * pli->Z;
        let evalue = first_domain.lnP.exp() * unsafe { (*hmmsearch_result.c_pli).Z };
        println!("First domain evalue: {:?}", evalue);

        assert_eq!(evalue, 1.4970530541655288e-48);
        assert_eq!(hmmsearch_result.nreported(), 1);
    }
}
