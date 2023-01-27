mod libhmmer_sys_extras;
mod hmm;
mod hmmsearch;

use log::*;
use std::ffi::CStr;
use std::fmt::Debug;

pub use crate::hmm::*;
pub use crate::hmmsearch::*;

pub enum Alphabet {
    Protein,
    RNA,
    DNA,
}

pub struct EaselSequence {
    // TODO: Implement Drop trait to free this
    pub c_sq: *mut libhmmer_sys::ESL_SQ,
}

impl EaselSequence {
    pub fn new(alphabet: Alphabet) -> Self {
        let c_alphabet = match alphabet {
            // *const libhmmer_sys::ESL_ALPHABET
            Alphabet::Protein => unsafe {
                libhmmer_sys::esl_alphabet_Create(libhmmer_sys::eslAMINO.try_into().unwrap())
            },
            Alphabet::RNA => unsafe {
                libhmmer_sys::esl_alphabet_Create(libhmmer_sys::eslRNA.try_into().unwrap())
            },
            Alphabet::DNA => unsafe {
                libhmmer_sys::esl_alphabet_Create(libhmmer_sys::eslDNA.try_into().unwrap())
            },
        };
        let c_sq = unsafe { libhmmer_sys::esl_sq_CreateDigital(c_alphabet) };
        Self { c_sq }
    }

    /// Replace (or initialise) the sequence data in this object with the given
    /// sequence. This method is a simplification of the sqascii_ReadSequence()
    /// C function in easel, in esl_sqio_ascii.c.
    ///
    /// The input sequence is assumed to be in the same alphabet as the one used
    /// to instantiate this struct. It is not NULL terminated.
    ///
    /// The seq given here is converted into a newly allocated dsq, so does not
    /// need to live after this function returns.
    pub fn replace_sequence(&mut self, seq: &[u8]) -> Result<(), &'static str> {
        let n = seq.len() as i64;

        unsafe {
            // free() previous sequence
            if (*self.c_sq).dsq.is_null() {
                debug!("Freeing previous dsq pointer");
                libc::free((*self.c_sq).dsq as *mut libc::c_void);
            }
            (*self.c_sq).dsq = libc::malloc(seq.len() + 2) as *mut u8;

            // esl_abc_Digitize(const ESL_ALPHABET *a, const char *seq, ESL_DSQ *dsq)
            self.digitise_sequence(seq)?;

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
        Ok(())
    }

    // Reimplementation of libhmmer_sys::esl_abc_Digitize but don't require a
    // NULL terminated sequence as input. Assumes self.dsq is already allocated.
    fn digitise_sequence(&mut self, seq: &[u8]) -> Result<(), &'static str> {
        // let sstatus = libhmmer_sys::esl_abc_Digitize((*self.c_sq).abc, seq.as_ptr() as *const i8, (*self.c_sq).dsq);

        // int     status;
        // int64_t i;			/* position in seq */
        // int64_t j;			/* position in dsq */
        // ESL_DSQ x;

        // status = eslOK;
        // dsq[0] = eslDSQ_SENTINEL;
        // for (i = 0, j = 1; seq[i] != '\0'; i++)
        //   {
        //     x = a->inmap[(int) seq[i]];
        //     if      (esl_abc_XIsValid(a, x)) dsq[j] = x;
        //     else if (x == eslDSQ_IGNORED) continue;
        //     else {
        //   status   = eslEINVAL;
        //   dsq[j] = esl_abc_XGetUnknown(a);
        //     }
        //     j++;
        //   }
        // dsq[j] = eslDSQ_SENTINEL;
        // return status;

        // easel/esl_alphabet.h:#define esl_abc_XIsValid(a, x)       ((x) < (a)->Kp)

        // Set initial sentinal
        unsafe {
            *(*self.c_sq).dsq = libhmmer_sys::eslDSQ_SENTINEL as u8;
        };

        // Set actual sequence
        #[allow(non_snake_case)]
        let Kp: u8 = unsafe { (*(*self.c_sq).abc).Kp.try_into().unwrap() };
        for (i, s) in seq.iter().enumerate() {
            let x = unsafe { (*(*self.c_sq).abc).inmap[*s as usize] };
            if x < Kp {
                // easel/esl_alphabet.h:#define esl_abc_XGetUnknown(a)       ((a)->Kp)
                unsafe {
                    *(*self.c_sq).dsq.add(i + 1) = x;
                };
            } else if x == libhmmer_sys::eslDSQ_IGNORED as u8 {
                continue;
            } else {
                return Err("Invalid character in sequence");
            }
        }

        // Set final sentinal
        unsafe {
            *((*self.c_sq).dsq.offset(seq.len() as isize + 1)) =
                libhmmer_sys::eslDSQ_SENTINEL as u8;
        };

        Ok(())
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
                .field(
                    "dsq",
                    &CStr::from_ptr((*self.c_sq).dsq as *mut i8).to_string_lossy(),
                )
                .field("dsq ptr", &(*self.c_sq).dsq)
                .field(
                    "dsq length",
                    &libhmmer_sys::esl_abc_dsqlen((*self.c_sq).dsq),
                )
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

impl Drop for EaselSequence {
    fn drop(&mut self) {
        unsafe {
            libhmmer_sys::esl_sq_Destroy(self.c_sq);
        }
    }
}
