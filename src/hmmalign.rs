use std::ffi::CStr;

use crate::{hmm::Hmm, EaselSequence};

#[derive(Debug)]
pub enum HmmerAlignError {
    AlignmentFailure,
    SerializationFailure,
    NoSequences,
}

impl std::fmt::Display for HmmerAlignError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HmmerAlignError::AlignmentFailure => write!(f, "Alignment/trace computation failed"),
            HmmerAlignError::SerializationFailure => write!(f, "Stockholm serialization failed"),
            HmmerAlignError::NoSequences => write!(f, "No sequences provided"),
        }
    }
}

impl std::error::Error for HmmerAlignError {}

/// Owned wrapper around an ESL_MSA (multiple sequence alignment).
pub struct EaselMsa {
    c_msa: *mut libhmmer_sys::ESL_MSA,
}

impl EaselMsa {
    /// Number of sequences in the alignment.
    pub fn num_sequences(&self) -> usize {
        unsafe { (*self.c_msa).nseq as usize }
    }

    /// Alignment length (number of columns, including gap columns).
    pub fn alignment_length(&self) -> usize {
        unsafe { (*self.c_msa).alen as usize }
    }

    /// Get the name of the i-th sequence.
    pub fn sequence_name(&self, i: usize) -> String {
        assert!(i < self.num_sequences(), "sequence index out of bounds");
        unsafe {
            let name_ptr = *(*self.c_msa).sqname.add(i);
            CStr::from_ptr(name_ptr).to_string_lossy().into_owned()
        }
    }

    /// Get the aligned sequence string for the i-th sequence (text mode).
    /// Returns the alignment row including gap characters.
    pub fn aligned_sequence(&self, i: usize) -> String {
        assert!(i < self.num_sequences(), "sequence index out of bounds");
        unsafe {
            // Check if we have text-mode alignment (aseq) or digital (ax)
            if !(*self.c_msa).aseq.is_null() {
                let seq_ptr = *(*self.c_msa).aseq.add(i);
                CStr::from_ptr(seq_ptr).to_string_lossy().into_owned()
            } else {
                // Digital mode — need to textize first, or decode manually
                // p7_tracealign_Seqs produces text-mode MSAs, so aseq should be set
                panic!("MSA is in digital mode; text mode expected");
            }
        }
    }

    /// Serialize the MSA to Stockholm format.
    pub fn to_stockholm(&self) -> Result<String, HmmerAlignError> {
        unsafe { msa_to_stockholm(self.c_msa) }
    }
}

impl Drop for EaselMsa {
    fn drop(&mut self) {
        unsafe {
            libhmmer_sys::esl_msa_Destroy(self.c_msa);
        }
    }
}

pub struct HmmerAlign {
    hmm: *mut libhmmer_sys::P7_HMM,
}

impl HmmerAlign {
    /// Create an aligner for the given HMM.
    pub fn new(hmm: &Hmm) -> Self {
        HmmerAlign { hmm: hmm.c_hmm }
    }

    /// Align a slice of sequences to the HMM and return the MSA as an
    /// owned [`EaselMsa`].
    ///
    /// Uses the same algorithm as the `hmmalign` command-line tool:
    /// optimal accuracy alignment via Forward/Backward/Decoding.
    pub fn align_sequences(
        &self,
        sequences: &[EaselSequence],
    ) -> Result<EaselMsa, HmmerAlignError> {
        let nseq = sequences.len();
        if nseq == 0 {
            return Err(HmmerAlignError::NoSequences);
        }

        unsafe {
            let hmm_m = (*self.hmm).M;

            let mut sq_ptrs: Vec<*mut libhmmer_sys::ESL_SQ> =
                sequences.iter().map(|s| s.c_sq).collect();

            let mut traces: Vec<*mut libhmmer_sys::P7_TRACE> = (0..nseq)
                .map(|_| libhmmer_sys::p7_trace_CreateWithPP())
                .collect();

            let trace_status = libhmmer_sys::p7_tracealign_computeTraces(
                self.hmm,
                sq_ptrs.as_mut_ptr(),
                0,
                nseq as i32,
                traces.as_mut_ptr(),
            );
            if trace_status != libhmmer_sys::eslOK as i32 {
                for t in &traces {
                    libhmmer_sys::p7_trace_Destroy(*t);
                }
                return Err(HmmerAlignError::AlignmentFailure);
            }

            let mut msa: *mut libhmmer_sys::ESL_MSA = std::ptr::null_mut();
            let align_status = libhmmer_sys::p7_tracealign_Seqs(
                sq_ptrs.as_mut_ptr(),
                traces.as_mut_ptr(),
                nseq as i32,
                hmm_m,
                crate::libhmmer_sys_extras::p7_DEFAULT,
                self.hmm,
                &mut msa,
            );

            for t in &traces {
                libhmmer_sys::p7_trace_Destroy(*t);
            }

            if align_status != libhmmer_sys::eslOK as i32 {
                return Err(HmmerAlignError::AlignmentFailure);
            }

            Ok(EaselMsa { c_msa: msa })
        }
    }

    /// Align a slice of sequences to the HMM and return the MSA
    /// in Stockholm format as a String.
    pub fn align_sequences_into_stockholm(
        &self,
        sequences: &[EaselSequence],
    ) -> Result<String, HmmerAlignError> {
        let msa = self.align_sequences(sequences)?;
        msa.to_stockholm()
    }
}

/// Serialize an ESL_MSA to Stockholm format string using open_memstream.
unsafe fn msa_to_stockholm(msa: *mut libhmmer_sys::ESL_MSA) -> Result<String, HmmerAlignError> {
    let mut buf: *mut libc::c_char = std::ptr::null_mut();
    let mut buf_size: libc::size_t = 0;
    let fp = libc::open_memstream(&mut buf, &mut buf_size);
    if fp.is_null() {
        return Err(HmmerAlignError::SerializationFailure);
    }

    let write_status = libhmmer_sys::esl_msafile_Write(
        fp as *mut libhmmer_sys::FILE,
        msa,
        libhmmer_sys::eslMSAFILE_STOCKHOLM as i32,
    );
    libc::fclose(fp);

    if write_status != libhmmer_sys::eslOK as i32 {
        if !buf.is_null() {
            libc::free(buf as *mut libc::c_void);
        }
        return Err(HmmerAlignError::SerializationFailure);
    }

    let result = if !buf.is_null() && buf_size > 0 {
        CStr::from_ptr(buf).to_string_lossy().into_owned()
    } else {
        String::new()
    };

    if !buf.is_null() {
        libc::free(buf as *mut libc::c_void);
    }

    Ok(result)
}
