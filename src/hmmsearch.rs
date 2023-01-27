use log::*;
use std::ffi::{CStr, CString};

use crate::{hmm::*, libhmmer_sys_extras::*, EaselSequence};

pub struct HmmerPipeline {
    info: HmmsearchWorkerInfo,
}

impl HmmerPipeline {
    pub fn new(hmm: &Hmm) -> HmmerPipeline {
        // P7_PROFILE      *gm      = NULL;
        // P7_OPROFILE     *om      = NULL;       /* optimized query profile                  */
        #[allow(unused_assignments)]
        let mut gm: *mut libhmmer_sys::P7_PROFILE = std::ptr::null_mut();
        #[allow(unused_assignments)]
        let mut om: *mut libhmmer_sys::P7_OPROFILE = std::ptr::null_mut();
        // int              textw    = 0;

        // WORKER_INFO     *info     = NULL;
        let bg = unsafe { libhmmer_sys::p7_bg_Create(hmm.c_alphabet()) };
        debug!("Background model created successfully");

        //   /* Convert to an optimized model */
        // gm = p7_profile_Create (hmm->M, abc);
        // om = p7_oprofile_Create(hmm->M, abc);
        // p7_ProfileConfig(hmm, info->bg, gm, 100, p7_LOCAL); /* 100 is a dummy length for now; and MSVFilter requires local mode */
        // p7_oprofile_Convert(gm, om);                  /* <om> is now p7_LOCAL, multihit */
        let abc = hmm.c_alphabet();
        gm = unsafe { libhmmer_sys::p7_profile_Create(hmm.length() as i32, abc) };
        debug!("Profile created successfully");
        om = unsafe { libhmmer_sys::p7_oprofile_Create(hmm.length() as i32, abc) };
        debug!("Optimized profile created successfully");

        unsafe {
            libhmmer_sys::p7_ProfileConfig(hmm.c_hmm, bg, gm, 100, p7_LOCAL);
            debug!("Profile configured successfully");
            libhmmer_sys::p7_oprofile_Convert(gm, om);
            debug!("Optimized profile converted successfully");
        }

        // /* Create processing pipeline and hit list */
        // info[i].th  = p7_tophits_Create();
        // info[i].om  = p7_oprofile_Clone(om);
        // info[i].pli = p7_pipeline_Create(go, om->M, 100, FALSE, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
        // status = p7_pli_NewModel(info[i].pli, info[i].om, info[i].bg);
        // if (status == eslEINVAL) p7_Fail(info->pli->errbuf);
        let th = unsafe { libhmmer_sys::p7_tophits_Create() };
        debug!("Tophits created successfully");
        // let om = libhmmer_sys::p7_oprofile_Clone(om);
        let pli = unsafe {
            libhmmer_sys::p7_pipeline_Create(
                std::ptr::null_mut(),
                (*om).M,
                100,
                0,
                p7_SEARCH_SEQS as u32,
            )
        };
        debug!("Pipeline created successfully");
        unsafe {
            let status = libhmmer_sys::p7_pli_NewModel(pli, om, bg);
            if status == eslEINVAL {
                panic!(); // TODO: Better msg
            }
        }
        debug!("Pipeline new model created successfully");

        let info = HmmsearchWorkerInfo {
            th: th,
            om: om,
            pli: pli,
            bg: bg,
        };

        HmmerPipeline { info: info }
    }

    pub fn run_hmm_on_file(&mut self, hmm: &Hmm, fasta_path: &std::path::Path) -> HmmsearchResult {
        debug!("Starting run_hmm_on_file");
        #[allow(unused_mut)]
        let mut dbfile = Self::open_target_sequences(&fasta_path.to_string_lossy());
        debug!("Target sequences opened successfully");

        // TODO: output_header(ofp, go, cfg->hmmfile, cfg->dbfile);

        //       esl_sqfile_SetDigital(dbfp, abc); //ReadBlock requires knowledge of the alphabet to decide how best to read blocks
        unsafe {
            libhmmer_sys::esl_sqfile_SetDigital(dbfile, hmm.c_alphabet());
        }
        debug!("Target sequences set to digital successfully");

        // if (fprintf(ofp, "Query:       %s  [M=%d]\n", hmm->name, hmm->M)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
        // if (hmm->acc)  { if (fprintf(ofp, "Accession:   %s\n", hmm->acc)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
        // if (hmm->desc) { if (fprintf(ofp, "Description: %s\n", hmm->desc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
        println!("Query:       {}  [M={}]", hmm.name(), hmm.length());
        println!("Accession:   {}", hmm.acc());
        println!("Description: {}", hmm.desc());

        // sstatus = serial_loop(info, dbfp, cfg->n_targetseq);
        // TODO: n_targetseq == -1 means no limit. OK for now.
        let sstatus = self.serial_loop_over_esl_sqio(dbfile, -1);

        // switch(sstatus)
        // {
        // case eslEFORMAT:
        //   esl_fatal("Parse failed (sequence file %s):\n%s\n",
        //       dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
        //   break;
        // case eslEOF:
        //   /* do nothing */
        //   break;
        // default:
        //   esl_fatal("Unexpected error %d reading sequence file %s", sstatus, dbfp->filename);
        // }
        #[allow(non_upper_case_globals)]
        match sstatus {
            eslEFORMAT => {
                // panic!("Parse failed (sequence file {}):\n{}", fasta_path.to_string_lossy(), unsafe { libhmmer_sys::esl_sqfile_GetErrorBuf(dbfile) });
                // TODO: Make the above compile
                panic!(
                    "Parse failed (sequence file {})",
                    fasta_path.to_string_lossy()
                );
            }
            eslEOF => {
                // do nothing
            }
            _ => {
                panic!(
                    "Unexpected error {} reading sequence file {}",
                    sstatus,
                    fasta_path.to_string_lossy()
                );
            }
        }

        // /* merge the results of the search results */
        // for (i = 1; i < infocnt; ++i)
        // {
        // p7_tophits_Merge(info[0].th, info[i].th);
        // p7_pipeline_Merge(info[0].pli, info[i].pli);

        // p7_pipeline_Destroy(info[i].pli);
        // p7_tophits_Destroy(info[i].th);
        // p7_oprofile_Destroy(info[i].om);
        // }

        // ===> This block of code is not necessary because we only have one thread.

        // /* Print the results.  */
        // p7_tophits_SortBySortkey(info->th);
        // p7_tophits_Threshold(info->th, info->pli);
        // p7_tophits_Targets(ofp, info->th, info->pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
        // p7_tophits_Domains(ofp, info->th, info->pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
        unsafe {
            libhmmer_sys::p7_tophits_SortBySortkey(self.info.th);
            libhmmer_sys::p7_tophits_Threshold(self.info.th, self.info.pli);

            // TODO: Implement (optional) output of default output.
            // panic!("Need to get FILE* from fopen or stdout to do this");
            // libhmmer_sys::p7_tophits_Targets(FIXME, (*info).th, (*info).pli, textw);
            // libhmmer_sys::p7_tophits_Domains(FIXME, (*info).th, (*info).pli, textw);
        }

        // if (tblfp)     p7_tophits_TabularTargets(tblfp,    hmm->name, hmm->acc, info->th, info->pli, (nquery == 1));
        // if (domtblfp)  p7_tophits_TabularDomains(domtblfp, hmm->name, hmm->acc, info->th, info->pli, (nquery == 1));
        // if (pfamtblfp) p7_tophits_TabularXfam(pfamtblfp, hmm->name, hmm->acc, info->th, info->pli);
        // TODO: The above. I don't need them for now.

        // TODO: Destroy, free, etc.

        return HmmsearchResult {
            c_th: self.info.th,
            c_pli: self.info.pli,
        };
    }

    // TODO: Add a coresponding free() method
    fn open_target_sequences(fasta_file: &str) -> *mut libhmmer_sys::esl_sqio_s {
        //   int              dbfmt    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */
        let dbfmt = 0; //libhmmer_sys::eslSQFILE_UNKNOWN;

        //   ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */
        let mut dbfp: *mut libhmmer_sys::ESL_SQFILE = std::ptr::null_mut();

        //   /* Open the target sequence database */
        //   status = esl_sqfile_Open(cfg->dbfile, dbfmt, p7_SEQDBENV, &dbfp);
        //   if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",          cfg->dbfile);
        //   else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",            cfg->dbfile);
        //   else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
        //   else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, cfg->dbfile);
        let file_pointer = CString::new(fasta_file.as_bytes()).unwrap().into_raw();
        let status = unsafe {
            // Open the file not assuming anything about its format, and let
            // autodetect do its thing. Possibly we should use eslSQFILE_FASTA.
            libhmmer_sys::esl_sqfile_Open(file_pointer, dbfmt, std::ptr::null(), &mut dbfp)
        };
        println!("Opened fasta file with status {}", status);

        if status == eslENOTFOUND {
            panic!("Failed to open sequence file {} for reading", fasta_file);
        } else if status == eslEFORMAT {
            panic!("Sequence file {} is empty or misformatted", fasta_file);
        } else if status == eslEINVAL {
            panic!("Can't autodetect format of a stdin or .gz seqfile");
        } else if status != eslOK {
            panic!(
                "Unexpected error {} opening sequence file {}",
                status, fasta_file
            );
        }
        return dbfp;
    }

    /// This method (called serial_loop in C) is not available in libhmmer_sys,
    /// so we have to implement it here. Intended as a direct replacement for
    /// the C function.
    fn serial_loop_over_esl_sqio(
        &mut self,
        dbfp: *mut libhmmer_sys::esl_sqio_s,
        n_targetseqs: i32,
    ) -> i32 {
        debug!("serial_loop");

        //   int              status;                       /* easel return code                               */
        let mut sstatus: i32;

        //   ESL_SQ   *dbsq     = NULL;   /* one target sequence (digital)  */
        #[allow(unused_mut, unused_assignments)]
        let mut dbsq: *mut libhmmer_sys::ESL_SQ = std::ptr::null_mut();

        // int seq_cnt = 0;
        let mut seq_cnt: i32 = 0;

        // For convenience
        let info = &mut self.info;

        // dbsq = esl_sq_CreateDigital(info->om->abc);
        dbsq = unsafe {
            let abc_here = (*info.om).abc;
            debug!("Creating digital sequence from abc {:?}", abc_here);
            debug!("K in abc is {}", (*abc_here).K);
            libhmmer_sys::esl_sq_CreateDigital(abc_here)
        };
        debug!("dbsq created as digital, from abc {:?}", unsafe {
            (*dbsq).abc
        });

        //   /* Main loop: */
        //   while ( (n_targetseqs==-1 || seq_cnt<n_targetseqs) &&  (sstatus = esl_sqio_Read(dbfp, dbsq)) == eslOK)
        //   {
        sstatus = unsafe { libhmmer_sys::esl_sqio_Read(dbfp, dbsq) };
        debug!("esl_sqio_Read returned {}", sstatus);
        debug!("dbsq is {:?}", dbsq);
        debug!("dbsq internals: {:#?}", EaselSequence { c_sq: dbsq });

        while (n_targetseqs == -1 || seq_cnt < n_targetseqs) && sstatus == eslOK {
            unsafe {
                // p7_pli_NewSeq(info->pli, dbsq);
                libhmmer_sys::p7_pli_NewSeq(info.pli, dbsq);

                // p7_bg_SetLength(info->bg, dbsq->n);
                libhmmer_sys::p7_bg_SetLength(
                    info.bg,
                    (*dbsq).n.try_into().expect("i64 -> i32 failed"),
                );
                // p7_oprofile_ReconfigLength(info->om, dbsq->n);
                libhmmer_sys::p7_oprofile_ReconfigLength(
                    info.om,
                    (*dbsq).n.try_into().expect("i64 -> i32 failed"),
                );

                // p7_Pipeline(info->pli, info->om, info->bg, dbsq, NULL, info->th);
                let p7_sstatus = libhmmer_sys::p7_Pipeline(
                    info.pli,
                    info.om,
                    info.bg,
                    dbsq,
                    std::ptr::null_mut(),
                    info.th,
                );
                if p7_sstatus != eslOK {
                    panic!("p7_Pipeline sstatus indicated failure, was {}", p7_sstatus);
                }

                // In the C code, this is part of the while loop condition.
                sstatus = libhmmer_sys::esl_sqio_Read(dbfp, dbsq);
                debug!("esl_sqio_Read returned {}", sstatus);
            }

            // seq_cnt++;
            seq_cnt += 1;
            // esl_sq_Reuse(dbsq);
            // p7_pipeline_Reuse(info->pli);
            unsafe {
                libhmmer_sys::esl_sq_Reuse(dbsq);
                libhmmer_sys::p7_pipeline_Reuse(info.pli);
            }
        }

        // if (n_targetseqs!=-1 && seq_cnt==n_targetseqs)
        // sstatus = eslEOF;
        if n_targetseqs != -1 && seq_cnt == n_targetseqs {
            sstatus = eslEOF;
        }

        // esl_sq_Destroy(dbsq);
        unsafe {
            libhmmer_sys::esl_sq_Destroy(dbsq);
        }

        return sstatus;
    }

    pub fn query(&mut self, easel_sequence: &EaselSequence) {
        let info = &mut self.info;

        unsafe {
            // p7_pli_NewSeq(info->pli, dbsq);
            if libhmmer_sys::p7_pli_NewSeq(info.pli, easel_sequence.c_sq) != eslOK {
                panic!()
            };

            // p7_bg_SetLength(info->bg, dbsq->n);
            if libhmmer_sys::p7_bg_SetLength(
                info.bg,
                (*easel_sequence.c_sq)
                    .n
                    .try_into()
                    .expect("i64 -> i32 failed"),
            ) != eslOK
            {
                panic!()
            };
            // p7_oprofile_ReconfigLength(info->om, dbsq->n);
            if libhmmer_sys::p7_oprofile_ReconfigLength(
                info.om,
                (*easel_sequence.c_sq)
                    .n
                    .try_into()
                    .expect("i64 -> i32 failed"),
            ) != eslOK
            {
                panic!()
            };

            // p7_Pipeline(info->pli, info->om, info->bg, dbsq, NULL, info->th);
            let sstatus = libhmmer_sys::p7_Pipeline(
                info.pli,
                info.om,
                info.bg,
                easel_sequence.c_sq,
                std::ptr::null_mut(),
                info.th,
            );
            debug!("query p7_Pipeline sstatus {}", sstatus);
            if sstatus != eslOK {
                panic!("p7_Pipeline sstatus indicated failure, was {}", sstatus);
            }
        }
    }

    pub fn get_results(&mut self) -> HmmsearchResult {
        unsafe {
            debug!("Running p7_tophits_SortBySortkey");
            libhmmer_sys::p7_tophits_SortBySortkey(self.info.th);
            debug!("Running p7_tophits_Threshold");
            libhmmer_sys::p7_tophits_Threshold(self.info.th, self.info.pli);
        }
        HmmsearchResult {
            c_th: self.info.th,
            c_pli: self.info.pli,
        }
    }
}

// typedef struct {
//       P7_BG            *bg;	         /* null model                              */
//       P7_PIPELINE      *pli;         /* work pipeline                           */
//       P7_TOPHITS       *th;          /* top hit results                         */
//       P7_OPROFILE      *om;          /* optimized query profile                 */
//     } WORKER_INFO;
pub struct HmmsearchWorkerInfo {
    bg: *mut libhmmer_sys::P7_BG,
    pli: *mut libhmmer_sys::P7_PIPELINE,
    th: *mut libhmmer_sys::P7_TOPHITS,
    om: *mut libhmmer_sys::P7_OPROFILE,
}

#[derive(Debug)]
pub struct HmmsearchResult {
    pub c_th: *mut libhmmer_sys::p7_tophits_s,
    pub c_pli: *mut libhmmer_sys::p7_pipeline_s,
}

impl Drop for HmmsearchResult {
    fn drop(&mut self) {
        unsafe {
            libhmmer_sys::p7_tophits_Destroy(self.c_th);
            libhmmer_sys::p7_pipeline_Destroy(self.c_pli);
        }
    }
}

impl HmmsearchResult {
    /// Number of reported hits
    pub fn nreported(&self) -> usize {
        unsafe {
            (*self.c_th)
                .nreported
                .try_into()
                .expect("i64 -> usize failed")
        }
    }

    pub fn hits(&self) -> HmmsearchResultTopHits {
        HmmsearchResultTopHits {
            c_th: self.c_th,
            c_pli: self.c_pli,
            current: 0,
            nreported: self.nreported(),
        }
    }
}

pub struct HmmsearchResultTopHits {
    c_th: *mut libhmmer_sys::p7_tophits_s,
    c_pli: *mut libhmmer_sys::p7_pipeline_s,
    current: usize,
    nreported: usize,
}

impl Iterator for HmmsearchResultTopHits {
    type Item = HmmsearchResultTopHit;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current >= self.nreported {
            return None;
        }

        let hit = HmmsearchResultTopHit {
            c_hit: unsafe { *((*self.c_th).hit.offset(self.current as isize)) },
            c_pli: self.c_pli,
            current_domain: 0,
        };

        self.current += 1;

        Some(hit)
    }
}

pub struct HmmsearchResultTopHit {
    c_hit: *mut libhmmer_sys::p7_hit_s,
    c_pli: *mut libhmmer_sys::p7_pipeline_s,
    current_domain: usize,
}

impl HmmsearchResultTopHit {
    // println!("Name of first hit {}", unsafe {
    //     CStr::from_ptr(first_hit.name).to_string_lossy()
    // });
    pub fn name(&self) -> String {
        unsafe {
            CStr::from_ptr((*self.c_hit).name)
                .to_string_lossy()
                .into_owned()
        }
    }

    // println!("Score of first hit overall {}", first_hit.score);
    pub fn score(&self) -> f32 {
        unsafe { (*self.c_hit).score }
    }
}

impl Iterator for HmmsearchResultTopHit {
    type Item = HmmsearchResultDomain;

    fn next(&mut self) -> Option<Self::Item> {
        // TODO: Check for end of hits
        if self.current_domain >= unsafe { (*self.c_hit).ndom.try_into().unwrap() } {
            return None;
        }
        println!("current domain counter {} ", self.current_domain);
        let domain = HmmsearchResultDomain {
            c_dom: unsafe { (*self.c_hit).dcl.offset(self.current_domain as isize) },
            c_pli: self.c_pli,
        };

        self.current_domain += 1;

        Some(domain)
    }
}

pub struct HmmsearchResultDomain {
    c_dom: *mut libhmmer_sys::p7_dom_s,
    c_pli: *mut libhmmer_sys::p7_pipeline_s,
}

impl HmmsearchResultDomain {
    // println!("First domain score: {}", first_domain.bitscore);
    pub fn bitscore(&self) -> f32 {
        unsafe { (*self.c_dom).bitscore }
    }

    // // exp(th->hit[h]->lnP) * pli->Z;
    // let evalue = first_domain.lnP.exp() * unsafe { (*hmmsearch_result.c_pli).Z };
    // println!("First domain evalue: {:?}", evalue);
    pub fn evalue(&self) -> f64 {
        unsafe { (*self.c_dom).lnP.exp() * (*self.c_pli).Z }
    }
}
