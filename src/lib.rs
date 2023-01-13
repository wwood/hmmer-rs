mod libhmmer_sys_extras;

mod hmmsearch;
mod hmm;

pub use crate::libhmmer_sys_extras::*;

pub use crate::hmm::*;
pub use crate::hmmsearch::*;

pub struct EaselSequence {
    c_sq: *mut libhmmer_sys::ESL_SQ,
}

impl EaselSequence {
    pub fn new(alphabet: *const libhmmer_sys::ESL_ALPHABET) -> Self {
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
    /// 
    /// TODO: Can we do even less here, if all we need is to satisfy the
    /// hmmsearch pipeline?
    pub fn replace_sequence(&mut self, seq: &[u8]) {
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