use log::*;
use env_logger;

use std::ffi::CString;
use std::ffi::CStr;

use libhmmer_sys;
use hmmer_rs::libhmmer_sys_extras;

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

// Static defines of things not available from libhmmer-sys e.g. #defines
// #define p7_LOCAL     1		/* multihit local:  "fs" mode   */
#[allow(non_upper_case_globals)]
static p7_LOCAL: i32 = 1;
// easel/easel.h:#define eslERRBUFSIZE 128
#[allow(non_upper_case_globals)]
static eslERRBUFSIZE: i32 = 128;

struct Hmm {
    pub c_hmm: *mut libhmmer_sys::P7_HMM,
}

impl Hmm {
    pub fn read_one_from_path(path: &std::path::Path) -> Result<Hmm, &'static str> {
        // char          errbuf[eslERRBUFSIZE];
        #[allow(unused_mut)]
        let mut errbuf = CString::new(vec![1; eslERRBUFSIZE as usize]).unwrap().into_raw();

        let hmmfile = CString::new(path.to_string_lossy().as_bytes()).unwrap().into_raw();
        let mut hfp: *mut libhmmer_sys::P7_HMMFILE = std::ptr::null_mut();

        // Read in a HMM from a file
        let status1 = unsafe {
            libhmmer_sys::p7_hmmfile_OpenE(hmmfile, std::ptr::null_mut(), &mut hfp, errbuf)
        };
        // eslOK = 0
        if status1 != 0 {
            error!("Error in initial reading of HMM file");
            return Err("Error in initial reading of HMM file");
        }
        debug!("HMM file opened successfully");

        // ESL_ALPHABET *abc     = NULL;	/* alphabet (set from the HMM file)*/
        // Set to NULL to not force alphabet
        let mut abc: *mut libhmmer_sys::ESL_ALPHABET = std::ptr::null_mut();
        // P7_HMM       *hmm     = NULL;
        let mut hmm: *mut libhmmer_sys::P7_HMM = std::ptr::null_mut();

        // Read the first HMM
        let status2 = unsafe {
            libhmmer_sys::p7_hmmfile_Read(hfp, &mut abc, &mut hmm)
        };
        if status2 != 0 {
            error!("Error in reading first HMM");
            return Err("Error in reading first HMM");
        }
        debug!("First HMM read successfully");

        // retake pointer to free memory
        unsafe {
            let _ = CString::from_raw(errbuf);
            let _ = CString::from_raw(hmmfile);
        }

        return Ok(Hmm {
            c_hmm: hmm,
        })
    }

    pub fn c_alphabet(&self) -> *const libhmmer_sys::ESL_ALPHABET {
        return unsafe {
            (*self.c_hmm).abc
        }
    }

    pub fn name(&self) -> String {
        return unsafe {
            CStr::from_ptr((*self.c_hmm).name).to_string_lossy().to_string()
        }
    }

    pub fn length(&self) -> u32 {
        return unsafe {
            (*self.c_hmm).M as u32
        }
    }

    pub fn acc(&self) -> String {
        let my_acc = unsafe {
            (*self.c_hmm).acc
        };
        if my_acc.is_null() {
            return "".to_string(); // Otherwise we get a segfault
        } else {
            return unsafe {
                CStr::from_ptr(my_acc).to_string_lossy().to_string()
            }
        }
    }

    pub fn desc(&self) -> String {
        let my_desc = unsafe {
            (*self.c_hmm).desc
        };
        if my_desc.is_null() {
            return "".to_string(); // Otherwise we get a segfault
        } else {
            return unsafe {
                CStr::from_ptr(my_desc).to_string_lossy().to_string()
            }
        }
    }
}

impl Drop for Hmm {
    fn drop(&mut self) {
        unsafe {
            libhmmer_sys::p7_hmm_Destroy(self.c_hmm);
        }
    }
}

struct HmmerPipeline {}

impl HmmerPipeline {
    pub fn run_hmm_on_file(&self, hmm: &Hmm, fasta_path: &std::path::Path) -> HmmsearchResult {
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

        // P7_PROFILE      *gm      = NULL;
        // P7_OPROFILE     *om      = NULL;       /* optimized query profile                  */
        #[allow(unused_assignments)]
        let mut gm: *mut libhmmer_sys::P7_PROFILE = std::ptr::null_mut();
        #[allow(unused_assignments)]
        let mut om: *mut libhmmer_sys::P7_OPROFILE = std::ptr::null_mut();
        // int              textw    = 0;

        // WORKER_INFO     *info     = NULL;
        let bg = unsafe {libhmmer_sys::p7_bg_Create(hmm.c_alphabet())};
        debug!("Background model created successfully");

        // if (fprintf(ofp, "Query:       %s  [M=%d]\n", hmm->name, hmm->M)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
        // if (hmm->acc)  { if (fprintf(ofp, "Accession:   %s\n", hmm->acc)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
        // if (hmm->desc) { if (fprintf(ofp, "Description: %s\n", hmm->desc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
        println!("Query:       {}  [M={}]", hmm.name(), hmm.length());
        println!("Accession:   {}", hmm.acc());
        println!("Description: {}", hmm.desc());
 
        
        //   /* Convert to an optimized model */
        // gm = p7_profile_Create (hmm->M, abc);
        // om = p7_oprofile_Create(hmm->M, abc);
        // p7_ProfileConfig(hmm, info->bg, gm, 100, p7_LOCAL); /* 100 is a dummy length for now; and MSVFilter requires local mode */
        // p7_oprofile_Convert(gm, om);                  /* <om> is now p7_LOCAL, multihit */
        let abc = hmm.c_alphabet();
        gm = unsafe {
            libhmmer_sys::p7_profile_Create(hmm.length() as i32, abc)
        };
        debug!("Profile created successfully");
        om = unsafe {
            libhmmer_sys::p7_oprofile_Create(hmm.length() as i32, abc)
        };
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
        let th = unsafe {libhmmer_sys::p7_tophits_Create()};
        debug!("Tophits created successfully");
        // let om = libhmmer_sys::p7_oprofile_Clone(om);
        let pli = unsafe {
            libhmmer_sys::p7_pipeline_Create(std::ptr::null_mut(), (*om).M, 100, 0, libhmmer_sys_extras::p7_SEARCH_SEQS as u32)
        };
        debug!("Pipeline created successfully");
        unsafe {
            let status = libhmmer_sys::p7_pli_NewModel(pli, om, bg);
            if status == libhmmer_sys_extras::eslEINVAL {
                panic!(); // TODO: Better msg
            }
        }
        debug!("Pipeline new model created successfully");

        // sstatus = serial_loop(info, dbfp, cfg->n_targetseq);
        // TODO: n_targetseq == -1 means no limit. OK for now.
        let mut info  = HmmsearchWorkerInfo {
            th: th,
            om: om,
            pli: pli,
            bg: bg,
        };
        let sstatus = self.serial_loop_over_esl_sqio(&mut info, dbfile, -1);

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
        match sstatus {
            libhmmer_sys_extras::eslEFORMAT => {
                // panic!("Parse failed (sequence file {}):\n{}", fasta_path.to_string_lossy(), unsafe { libhmmer_sys::esl_sqfile_GetErrorBuf(dbfile) });
                // TODO: Make the above compile
                panic!("Parse failed (sequence file {})", fasta_path.to_string_lossy());
            },
            libhmmer_sys_extras::eslEOF => {
                // do nothing
            },
            _ => {
                panic!("Unexpected error {} reading sequence file {}", sstatus, fasta_path.to_string_lossy());
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
            libhmmer_sys::p7_tophits_SortBySortkey(info.th);
            libhmmer_sys::p7_tophits_Threshold(info.th, info.pli);
            
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
            c_th: info.th,
            c_pli: info.pli,
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

        if status == libhmmer_sys_extras::eslENOTFOUND {
            panic!("Failed to open sequence file {} for reading", fasta_file);
        } else if status == libhmmer_sys_extras::eslEFORMAT {
            panic!("Sequence file {} is empty or misformatted", fasta_file);
        } else if status == libhmmer_sys_extras::eslEINVAL {
            panic!("Can't autodetect format of a stdin or .gz seqfile");
        } else if status != libhmmer_sys_extras::eslOK {
            panic!("Unexpected error {} opening sequence file {}", status, fasta_file);
        }
        return dbfp;
    }

    /// This method (called serial_loop in C) is not available in libhmmer_sys,
    /// so we have to implement it here. Intended as a direct replacement for
    /// the C function.
    fn serial_loop_over_esl_sqio(&self, info: &mut HmmsearchWorkerInfo, dbfp: *mut libhmmer_sys::esl_sqio_s, n_targetseqs: i32) -> i32 {
        debug!("serial_loop");

        //   int              status;                       /* easel return code                               */
        let mut sstatus: i32;

        //   ESL_SQ   *dbsq     = NULL;   /* one target sequence (digital)  */
        #[allow(unused_mut, unused_assignments)]
        let mut dbsq: *mut libhmmer_sys::ESL_SQ = std::ptr::null_mut();
        
        // int seq_cnt = 0;
        let mut seq_cnt: i32 = 0;

        // dbsq = esl_sq_CreateDigital(info->om->abc);
        dbsq = unsafe {
            let abc_here = (*info.om).abc;
            debug!("Creating digital sequence from abc {:?}", abc_here);
            debug!("K in abc is {}", (*abc_here).K);
            libhmmer_sys::esl_sq_CreateDigital(abc_here)
        };
        debug!("dbsq created as digital, from abc {:?}", unsafe { (*dbsq).abc });

        //   /* Main loop: */
        //   while ( (n_targetseqs==-1 || seq_cnt<n_targetseqs) &&  (sstatus = esl_sqio_Read(dbfp, dbsq)) == eslOK)
        //   {
        sstatus = unsafe { libhmmer_sys::esl_sqio_Read(dbfp, dbsq) };
        debug!("esl_sqio_Read returned {}", sstatus);
        while (n_targetseqs == -1 || seq_cnt < n_targetseqs) && sstatus == libhmmer_sys_extras::eslOK {
            unsafe {
                // p7_pli_NewSeq(info->pli, dbsq);
                libhmmer_sys::p7_pli_NewSeq(info.pli, dbsq);

                // p7_bg_SetLength(info->bg, dbsq->n);
                libhmmer_sys::p7_bg_SetLength(info.bg, (*dbsq).n.try_into().expect("i64 -> i32 failed"));
                // p7_oprofile_ReconfigLength(info->om, dbsq->n);
                libhmmer_sys::p7_oprofile_ReconfigLength(info.om, (*dbsq).n.try_into().expect("i64 -> i32 failed"));

                // p7_Pipeline(info->pli, info->om, info->bg, dbsq, NULL, info->th);
                libhmmer_sys::p7_Pipeline(info.pli, info.om, info.bg, dbsq, std::ptr::null_mut(), info.th);

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
            sstatus = libhmmer_sys_extras::eslEOF;
        }

        // esl_sq_Destroy(dbsq);
        unsafe {
            libhmmer_sys::esl_sq_Destroy(dbsq);
        }

        return sstatus;
    }

    pub fn query(easel_sequence: &EaselSequence, info: &mut HmmsearchWorkerInfo) {
        unsafe {
            // p7_pli_NewSeq(info->pli, dbsq);
            libhmmer_sys::p7_pli_NewSeq(info.pli, easel_sequence.c_sq);

            // p7_bg_SetLength(info->bg, dbsq->n);
            libhmmer_sys::p7_bg_SetLength(info.bg, (*easel_sequence.c_sq).n.try_into().expect("i64 -> i32 failed"));
            // p7_oprofile_ReconfigLength(info->om, dbsq->n);
            libhmmer_sys::p7_oprofile_ReconfigLength(info.om, (*easel_sequence.c_sq).n.try_into().expect("i64 -> i32 failed"));

            // p7_Pipeline(info->pli, info->om, info->bg, dbsq, NULL, info->th);
            libhmmer_sys::p7_Pipeline(info.pli, info.om, info.bg, easel_sequence.c_sq, std::ptr::null_mut(), info.th);
        }
    }
}


// typedef struct {
//       P7_BG            *bg;	         /* null model                              */
//       P7_PIPELINE      *pli;         /* work pipeline                           */
//       P7_TOPHITS       *th;          /* top hit results                         */
//       P7_OPROFILE      *om;          /* optimized query profile                 */
//     } WORKER_INFO;
struct HmmsearchWorkerInfo {
    bg: *mut libhmmer_sys::P7_BG,
    pli: *mut libhmmer_sys::P7_PIPELINE,
    th: *mut libhmmer_sys::P7_TOPHITS,
    om: *mut libhmmer_sys::P7_OPROFILE,
}

#[derive(Debug)]
struct HmmsearchResult {
    c_th: *mut libhmmer_sys::p7_tophits_s,
    c_pli: *mut libhmmer_sys::p7_pipeline_s,
}

impl Drop for HmmsearchResult {
    fn drop(&mut self) {
        unsafe {
            libhmmer_sys::p7_tophits_Destroy(self.c_th);
            libhmmer_sys::p7_pipeline_Destroy(self.c_pli);
        }
    }
}

struct EaselSequence {
    c_sq: *mut libhmmer_sys::ESL_SQ,
}

impl EaselSequence {
    fn new(alphabet: *const libhmmer_sys::ESL_ALPHABET) -> Self {
        let c_sq = unsafe { libhmmer_sys::esl_sq_CreateDigital(alphabet) };
        Self { c_sq }
    }

    /// Replace (or initialise) the sequence data in this object with the given
    /// sequence. This method is a simplification of the sqascii_ReadSequence()
    /// C function in easel.
    /// 
    /// The seq given here should remain valid for the lifetime of this object,
    /// or until this method is called again, because the internal C pointer
    /// will be incorrect, causing undefined behavior.
    fn replace_sequence(&mut self, seq: &[u8]) {
        let n = seq.len() as i64;

        unsafe {
            // addbuf(sqfp, sq, n);
            (*self.c_sq).seq = seq.as_ptr() as *mut i8;

            // sq->start = 1;
            (*self.c_sq).start = 1;
            // sq->end   = sq->n;
            (*self.c_sq).end = n;
            // sq->C     = 0;
            (*self.c_sq).C = 0;
            // sq->W     = sq->n;
            (*self.c_sq).W = n;
            // sq->L     = sq->n;
            (*self.c_sq).L = n;
        }
    }
}