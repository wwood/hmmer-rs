mod libhmmer_sys_extras;

mod hmmsearch;
mod hmm;

use std::ffi::CStr;
use std::fmt::Debug;
use log::*;
use libc;

pub use crate::libhmmer_sys_extras::*;

pub use crate::hmm::*;
pub use crate::hmmsearch::*;

pub struct EaselSequence {
    // TODO: Implement Drop trait to free this
    pub c_sq: *mut libhmmer_sys::ESL_SQ,
}

impl EaselSequence {
    pub fn new(alphabet: *const libhmmer_sys::ESL_ALPHABET) -> Self {
        let c_sq = unsafe { libhmmer_sys::esl_sq_CreateDigital(alphabet) };
        Self { c_sq }
    }

    /// Replace (or initialise) the sequence data in this object with the given
    /// sequence. This method is a simplification of the sqascii_ReadSequence()
    /// C function in easel, in esl_sqio_ascii.c.
    /// 
    /// The seq given here should remain valid for the lifetime of this object,
    /// or until this method is called again, because the internal C pointer
    /// will be incorrect, causing undefined behavior.
    /// 
    /// TODO: Can we do even less here, if all we need is to satisfy the
    /// hmmsearch pipeline?
    pub fn replace_sequence(&mut self, seq: &[u8]) -> Result<(), &'static str> {
        let n = seq.len() as i64 - 1; // minus one for NULL terminator

        unsafe {
            // free() previous sequence
            if (*self.c_sq).dsq != std::ptr::null_mut() {
                debug!("Freeing previous dsq pointer");
                libc::free((*self.c_sq).dsq as *mut libc::c_void);
            }
            (*self.c_sq).dsq = libc::malloc(seq.len()+2) as *mut u8;
            // esl_abc_Digitize(const ESL_ALPHABET *a, const char *seq, ESL_DSQ *dsq)
            // TODO: Check return value
            let sstatus = libhmmer_sys::esl_abc_Digitize((*self.c_sq).abc, seq.as_ptr() as *const i8, (*self.c_sq).dsq);
            debug!("esl_abc_Digitize returned {}", sstatus);
            if sstatus != 0 {
                return Err("esl_abc_Digitize returned non-zero exitstatus {}, potentially indicating an invalid character in the sequence")
            }

            (*self.c_sq).n = n;

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

        debug!("Replaced sequence, now have {:#?}", self);
        return Ok(());
    }
}

impl Debug for EaselSequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // name: *mut c_char
        // acc: *mut c_char
        // desc: *mut c_char
        // tax_id: i32
        // seq: *mut c_char
        // dsq: *mut ESL_DSQ
        // ss: *mut c_char
        // n: i64
        // start: i64
        // end: i64
        // C: i64
        // W: i64
        // L: i64
        // source: *mut c_char
        // nalloc: c_int
        // aalloc: c_int
        // dalloc: c_int
        // salloc: i64
        // srcalloc: c_int
        // idx: i64
        // roff: off_t
        // hoff: off_t
        // doff: off_t
        // eoff: off_t
        // xr_tag: *mut *mut c_char
        // xr: *mut *mut c_char
        // nxr: c_int
        // abc: *const ESL_ALPHABET
        unsafe {
            f.debug_struct("EaselSequence")
                .field("c_sq", &self.c_sq)
                .field("name", &CStr::from_ptr((*self.c_sq).name).to_string_lossy())
                .field("acc", &CStr::from_ptr((*self.c_sq).acc).to_string_lossy())
                .field("dsq", &CStr::from_ptr((*self.c_sq).dsq as *mut i8).to_string_lossy())
                .field("dsq ptr", &(*self.c_sq).dsq)
                .field("dsq length", &libhmmer_sys::esl_abc_dsqlen((*self.c_sq).dsq))
                .field("tax_id", &(*self.c_sq).tax_id)
                // .field("seq", &CStr::from_ptr((*self.c_sq).seq).to_string_lossy())
                // .field("ss", &CStr::from_ptr((*self.c_sq).ss).to_string_lossy())
                .field("n", &(*self.c_sq).n)
                .field("start", &(*self.c_sq).start)
                .field("end", &(*self.c_sq).end)
                .field("C", &(*self.c_sq).C)
                .field("W", &(*self.c_sq).W)
                .field("L", &(*self.c_sq).L)

                .finish()
        }
    }
}