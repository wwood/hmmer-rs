use std::ffi::CString;
use std::ffi::CStr;

use libhmmer_sys;

fn main() {
    println!("Hello, world!");

    let hmm = Hmm::read_one_from_path(std::path::Path::new("test/data/DNGNGWU00030_mingle_output_good_seqs.hmm"));

    println!("HMM name: {}", unsafe { CStr::from_ptr((*hmm.c_hmm).name).to_string_lossy() });
}

struct Hmm {
    pub c_hmm: *mut libhmmer_sys::P7_HMM,
}

impl Hmm {
    pub fn read_one_from_path(path: &std::path::Path) -> Hmm {
        // char          errbuf[eslERRBUFSIZE];
        // easel/easel.h:#define eslERRBUFSIZE 128
        #[allow(unused_mut)]
        let mut errbuf = CString::new(vec![1; 128]).unwrap().into_raw();

        let hmmfile = CString::new(path.to_string_lossy().as_bytes()).unwrap().into_raw();
        let mut hfp: *mut libhmmer_sys::P7_HMMFILE = std::ptr::null_mut();

        // Read in a HMM from a file
        let status1 = unsafe {
            libhmmer_sys::p7_hmmfile_OpenE(hmmfile, std::ptr::null_mut(), &mut hfp, errbuf)
        };
        // eslOK = 0
        if status1 != 0 {
            panic!("Error in initial reading of HMM file");
        }
        println!("HMM file opened successfully");

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
            panic!("Error in reading first HMM");
        }
        println!("First HMM read successfully");

        // retake pointer to free memory
        unsafe {
            let _ = CString::from_raw(errbuf);
            let _ = CString::from_raw(hmmfile);
        }

        return Hmm {
            c_hmm: hmm,
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