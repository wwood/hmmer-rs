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
        ))
        .unwrap();
        let hmm = &hmms[0];

        println!("HMM name: {}", unsafe {
            CStr::from_ptr((*hmm.c_hmm).name).to_string_lossy()
        });

        let mut hmmsearch = HmmerPipeline::new(hmm);

        let hmmsearch_result = hmmsearch.run_hmm_on_file(
            hmm,
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
        ))
        .unwrap();
        let hmm = &hmms[0];

        println!("HMM name: {}", hmm.name());

        let mut hmmsearch = HmmerPipeline::new(hmm);

        let mut query_seq = EaselSequence::new(Alphabet::Protein);
        let seq: &[u8] =
            b"MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ";

        query_seq.replace_sequence(seq).unwrap();
        debug!("Query seq replaced;");

        hmmsearch.query(&query_seq);

        let hmmsearch_result = hmmsearch.get_results();

        debug!("HMMsearch result: {:?}", hmmsearch_result);

        // th->hit[h]->score,
        assert_eq!(1, hmmsearch_result.nreported());

        let mut num_hits = 0;
        for hit in hmmsearch_result.hits() {
            num_hits += 1;
            assert_eq!(150.01991, hit.score());
            let domains = hit.collect::<Vec<_>>();

            assert_eq!(1, domains.len());
            let first_domain = &domains[0];

            assert_eq!(149.90887, first_domain.bitscore());

            assert_eq!(1.4970530541655288e-48, first_domain.evalue());
        }

        assert_eq!(1, num_hits);
    }

    #[test]
    fn test_set_name() {
        let mut seq = EaselSequence::new(Alphabet::Protein);
        seq.replace_sequence(b"MVYSGPNAPIEVGN").unwrap();
        seq.set_name("my_sequence").unwrap();

        let name = unsafe {
            CStr::from_ptr((*seq.c_sq).name)
                .to_string_lossy()
                .to_string()
        };
        assert_eq!(name, "my_sequence");
    }

    #[test]
    fn test_set_name_appears_in_hits() {
        let hmms = Hmm::read_hmms_from_path(std::path::Path::new(
            "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm",
        ))
        .unwrap();
        let hmm = &hmms[0];

        let mut pipeline = HmmerPipeline::new(hmm);

        let mut query_seq = EaselSequence::new(Alphabet::Protein);
        query_seq
            .replace_sequence(
                b"MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
            )
            .unwrap();
        query_seq.set_name("test_seq_name").unwrap();

        pipeline.query(&query_seq);
        let result = pipeline.get_results();

        assert_eq!(1, result.nreported());
        for hit in result.hits() {
            assert_eq!("test_seq_name", hit.name());
        }
    }

    #[test]
    fn test_search_sequences() {
        let hmms = Hmm::read_hmms_from_path(std::path::Path::new(
            "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm",
        ))
        .unwrap();
        let hmm = &hmms[0];

        let mut pipeline = HmmerPipeline::new(hmm);

        let mut seq1 = EaselSequence::new(Alphabet::Protein);
        seq1.replace_sequence(
            b"MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
        )
        .unwrap();
        seq1.set_name("seq1").unwrap();

        let mut seq2 = EaselSequence::new(Alphabet::Protein);
        seq2.replace_sequence(b"AAAAAAAAAAAAAAAAAAAA").unwrap();
        seq2.set_name("seq2_no_hit").unwrap();

        let sequences = vec![seq1, seq2];
        let result = pipeline.search_sequences(&sequences);

        assert_eq!(1, result.nreported());
        for hit in result.hits() {
            assert_eq!("seq1", hit.name());
            assert_eq!(150.01991, hit.score());
        }
    }

    #[test]
    fn test_search_sequences_can_be_called_multiple_times() {
        let hmms = Hmm::read_hmms_from_path(std::path::Path::new(
            "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm",
        ))
        .unwrap();
        let hmm = &hmms[0];

        let mut pipeline = HmmerPipeline::new(hmm);

        let mut seq1 = EaselSequence::new(Alphabet::Protein);
        seq1.replace_sequence(
            b"MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
        )
        .unwrap();
        seq1.set_name("seq1").unwrap();

        let sequences = vec![seq1];

        // First call
        let result1 = pipeline.search_sequences(&sequences);
        assert_eq!(1, result1.nreported());
        drop(result1);

        // Second call on same pipeline should also work
        let result2 = pipeline.search_sequences(&sequences);
        assert_eq!(1, result2.nreported());
    }

    #[test]
    fn test_evalue_threshold_filters_hits() {
        let hmms = Hmm::read_hmms_from_path(std::path::Path::new(
            "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm",
        ))
        .unwrap();
        let hmm = &hmms[0];

        // With a very strict E-value threshold, the hit should still pass
        // (it has a very strong e-value ~1.5e-48)
        let mut pipeline = HmmerPipeline::new(hmm).with_dom_evalue(1e-40);

        let mut seq = EaselSequence::new(Alphabet::Protein);
        seq.replace_sequence(
            b"MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
        )
        .unwrap();
        seq.set_name("good_hit").unwrap();

        let result = pipeline.search_sequences(&[seq]);
        assert_eq!(1, result.nreported());
    }

    #[test]
    fn test_evalue_threshold_excludes_weak_hits() {
        let hmms = Hmm::read_hmms_from_path(std::path::Path::new(
            "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm",
        ))
        .unwrap();
        let hmm = &hmms[0];

        // With an impossibly strict E-value, nothing should pass
        let mut pipeline = HmmerPipeline::new(hmm).with_seq_evalue(1e-200);

        let mut seq = EaselSequence::new(Alphabet::Protein);
        seq.replace_sequence(
            b"MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
        )
        .unwrap();
        seq.set_name("test").unwrap();

        let result = pipeline.search_sequences(&[seq]);
        assert_eq!(0, result.nreported());
    }

    #[test]
    fn test_hit_evalue() {
        let hmms = Hmm::read_hmms_from_path(std::path::Path::new(
            "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm",
        ))
        .unwrap();
        let hmm = &hmms[0];

        let mut pipeline = HmmerPipeline::new(hmm);

        let mut seq = EaselSequence::new(Alphabet::Protein);
        seq.replace_sequence(
            b"MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
        )
        .unwrap();
        seq.set_name("test").unwrap();

        let result = pipeline.search_sequences(&[seq]);
        assert_eq!(1, result.nreported());

        for hit in result.hits() {
            // Hit-level evalue should be very small
            assert!(hit.evalue() < 1e-40, "evalue was {}", hit.evalue());
            // Bitscore should match score
            assert_eq!(hit.bitscore(), hit.score());
            assert!(hit.bitscore() > 100.0);
        }
    }

    #[test]
    fn test_hmmalign() {
        let hmms = Hmm::read_hmms_from_path(std::path::Path::new(
            "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm",
        ))
        .unwrap();
        let hmm = &hmms[0];

        let aligner = HmmerAlign::new(hmm);

        let mut seq = EaselSequence::new(Alphabet::Protein);
        seq.replace_sequence(
            b"MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
        )
        .unwrap();
        seq.set_name("test_seq").unwrap();

        let msa = aligner.align_sequences(&[seq]).unwrap();

        assert_eq!(msa.num_sequences(), 1);
        assert_eq!(msa.sequence_name(0), "test_seq");
        assert!(msa.alignment_length() > 0);

        // The aligned sequence should contain the original residues (possibly with gaps)
        let aligned = msa.aligned_sequence(0);
        assert!(aligned.len() > 0);
        // Removing gaps and uppercasing should recover the original sequence
        // (HMMER uses lowercase for insert-state residues)
        let ungapped: String = aligned
            .chars()
            .filter(|c| *c != '-' && *c != '.')
            .flat_map(|c| c.to_uppercase())
            .collect();
        assert_eq!(
            ungapped,
            "MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ"
        );
    }

    #[test]
    fn test_hmmalign_multiple_sequences() {
        let hmms = Hmm::read_hmms_from_path(std::path::Path::new(
            "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm",
        ))
        .unwrap();
        let hmm = &hmms[0];

        let aligner = HmmerAlign::new(hmm);

        let mut seq1 = EaselSequence::new(Alphabet::Protein);
        seq1.replace_sequence(
            b"MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
        )
        .unwrap();
        seq1.set_name("seq_one").unwrap();

        let mut seq2 = EaselSequence::new(Alphabet::Protein);
        seq2.replace_sequence(
            b"MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
        )
        .unwrap();
        seq2.set_name("seq_two").unwrap();

        let msa = aligner.align_sequences(&[seq1, seq2]).unwrap();

        assert_eq!(msa.num_sequences(), 2);
        assert_eq!(msa.sequence_name(0), "seq_one");
        assert_eq!(msa.sequence_name(1), "seq_two");
        // Both sequences are identical, so aligned sequences should match
        assert_eq!(msa.aligned_sequence(0), msa.aligned_sequence(1));
    }

    #[test]
    fn test_hmmalign_into_stockholm() {
        let hmms = Hmm::read_hmms_from_path(std::path::Path::new(
            "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm",
        ))
        .unwrap();
        let hmm = &hmms[0];

        let aligner = HmmerAlign::new(hmm);

        let mut seq = EaselSequence::new(Alphabet::Protein);
        seq.replace_sequence(
            b"MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
        )
        .unwrap();
        seq.set_name("test_seq").unwrap();

        let stockholm = aligner.align_sequences_into_stockholm(&[seq]).unwrap();

        assert!(
            stockholm.contains("# STOCKHOLM 1.0"),
            "Expected Stockholm header, got: {}",
            &stockholm[..100.min(stockholm.len())]
        );
        assert!(stockholm.contains("//"), "Expected Stockholm terminator");
        assert!(
            stockholm.contains("test_seq"),
            "Expected sequence name in output"
        );
    }

    #[test]
    fn test_hmmalign_to_stockholm_via_msa() {
        let hmms = Hmm::read_hmms_from_path(std::path::Path::new(
            "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm",
        ))
        .unwrap();
        let hmm = &hmms[0];

        let aligner = HmmerAlign::new(hmm);

        let mut seq = EaselSequence::new(Alphabet::Protein);
        seq.replace_sequence(
            b"MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
        )
        .unwrap();
        seq.set_name("test_seq").unwrap();

        let msa = aligner.align_sequences(&[seq]).unwrap();
        let stockholm = msa.to_stockholm().unwrap();

        assert!(stockholm.contains("# STOCKHOLM 1.0"));
        assert!(stockholm.contains("test_seq"));
        assert!(stockholm.contains("//"));
    }

    #[test]
    fn test_hmmalign_empty_sequences() {
        let hmms = Hmm::read_hmms_from_path(std::path::Path::new(
            "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm",
        ))
        .unwrap();
        let hmm = &hmms[0];

        let aligner = HmmerAlign::new(hmm);
        let result = aligner.align_sequences(&[]);
        assert!(result.is_err());
    }
}
