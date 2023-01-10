use log::*;

use std::ffi::CString;
use std::ffi::CStr;

use libhmmer_sys;
use hmmer_rs::libhmmer_sys_extras;

fn main() {
    let hmm = Hmm::read_one_from_path(std::path::Path::new("test/data/DNGNGWU00030_mingle_output_good_seqs.hmm")).unwrap();

    println!("HMM name: {}", unsafe { CStr::from_ptr((*hmm.c_hmm).name).to_string_lossy() });

    HmmerPipeline::run_hmm_on_file(&hmm, std::path::Path::new("test/data/DNGNGWU00030_mingle_output_good_seqs.fasta"));
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
        return unsafe {
            CStr::from_ptr((*self.c_hmm).acc).to_string_lossy().to_string()
        }
    }

    pub fn desc(&self) -> String {
        return unsafe {
            CStr::from_ptr((*self.c_hmm).desc).to_string_lossy().to_string()
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
    pub fn run_hmm_on_file(hmm: &Hmm, fasta_path: &std::path::Path) {
        let mut dbfile = Self::open_target_sequences(&fasta_path.to_string_lossy());

        // P7_PROFILE      *gm      = NULL;
        // P7_OPROFILE     *om      = NULL;       /* optimized query profile                  */
        let mut gm: *mut libhmmer_sys::P7_PROFILE = std::ptr::null_mut();
        let mut om: *mut libhmmer_sys::P7_OPROFILE = std::ptr::null_mut();

        // WORKER_INFO     *info     = NULL;
        let mut info: *mut HmmsearchWorkerInfo = std::ptr::null_mut();

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
        om = unsafe {
            libhmmer_sys::p7_oprofile_Create(hmm.length() as i32, abc)
        };


        unsafe {
            // TODO: the null model here is wrong.
            libhmmer_sys::p7_ProfileConfig(hmm.c_hmm, std::ptr::null_mut(), gm, 100, p7_LOCAL);
            libhmmer_sys::p7_oprofile_Convert(gm, om);
            panic!();
        }
    }
    
    // TODO: Add a coresponding free() method
    fn open_target_sequences(fasta_file: &str) -> *mut libhmmer_sys::esl_sqio_s {
        //   int              dbfmt    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */
        let dbfmt = 0; //libhmmer_sys::eslSQFILE_UNKNOWN;

        // #define p7_SEQDBENV          "BLASTDB"
        #[allow(non_snake_case)]
        let p7_SEQDBENV = "BLASTDB";

        //   ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */
        let mut dbfp: *mut libhmmer_sys::ESL_SQFILE = std::ptr::null_mut();

        //   /* Open the target sequence database */
        //   status = esl_sqfile_Open(cfg->dbfile, dbfmt, p7_SEQDBENV, &dbfp);
        //   if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",          cfg->dbfile);
        //   else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",            cfg->dbfile);
        //   else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
        //   else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, cfg->dbfile); 
        let status = unsafe {
            libhmmer_sys::esl_sqfile_Open(fasta_file.as_ptr() as *const i8, dbfmt, p7_SEQDBENV.as_ptr() as *const i8, &mut dbfp)
        }; 

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