use pyo3::prelude::*;

use crate::{Alphabet, EaselSequence, Hmm, HmmerAlign, HmmerPipeline};

/// Convert a list of (name, sequence) Python tuples into EaselSequences.
fn tuples_to_easel_sequences(sequences: Vec<(String, String)>) -> PyResult<Vec<EaselSequence>> {
    let mut easel_seqs = Vec::with_capacity(sequences.len());
    for (name, seq) in &sequences {
        let mut es = EaselSequence::new(Alphabet::Protein);
        es.replace_sequence(seq.as_bytes())
            .map_err(pyo3::exceptions::PyValueError::new_err)?;
        es.set_name(name)
            .map_err(pyo3::exceptions::PyValueError::new_err)?;
        easel_seqs.push(es);
    }
    Ok(easel_seqs)
}

#[pyclass(name = "Hmm", unsendable)]
pub struct PyHmm {
    inner: Hmm,
}

#[pymethods]
impl PyHmm {
    /// Read all HMMs from a file, returning a list of Hmm objects.
    #[staticmethod]
    fn read_hmms_from_path(path: &str) -> PyResult<Vec<PyHmm>> {
        let hmms = Hmm::read_hmms_from_path(std::path::Path::new(path))
            .map_err(pyo3::exceptions::PyIOError::new_err)?;
        Ok(hmms.into_iter().map(|h| PyHmm { inner: h }).collect())
    }

    #[getter]
    fn name(&self) -> String {
        self.inner.name()
    }

    #[getter]
    fn length(&self) -> u32 {
        self.inner.length()
    }
}

#[pyclass(name = "HmmsearchHit")]
pub struct PyHmmsearchHit {
    #[pyo3(get)]
    name: String,
    #[pyo3(get)]
    score: f32,
    #[pyo3(get)]
    evalue: f64,
}

#[pyclass(name = "HmmerPipeline", unsendable)]
pub struct PyHmmerPipeline {
    inner: HmmerPipeline,
}

#[pymethods]
impl PyHmmerPipeline {
    #[new]
    fn new(hmm: &PyHmm) -> Self {
        PyHmmerPipeline {
            inner: HmmerPipeline::new(&hmm.inner),
        }
    }

    /// Set the per-sequence E-value threshold.
    fn with_seq_evalue(mut slf: PyRefMut<'_, Self>, e: f64) -> PyRefMut<'_, Self> {
        unsafe {
            (*slf.inner.pli()).E = e;
            (*slf.inner.pli()).use_bit_cutoffs = 0;
            (*slf.inner.pli()).inc_by_E = 1;
        }
        slf
    }

    /// Set the per-domain E-value threshold.
    fn with_dom_evalue(mut slf: PyRefMut<'_, Self>, e: f64) -> PyRefMut<'_, Self> {
        unsafe {
            (*slf.inner.pli()).domE = e;
            (*slf.inner.pli()).use_bit_cutoffs = 0;
            (*slf.inner.pli()).inc_by_E = 1;
        }
        slf
    }

    /// Run hmmsearch against a list of (name, sequence) tuples.
    fn search_sequences(
        &mut self,
        sequences: Vec<(String, String)>,
    ) -> PyResult<Vec<PyHmmsearchHit>> {
        let easel_seqs = tuples_to_easel_sequences(sequences)?;
        let result = self.inner.search_sequences(&easel_seqs);
        let mut hits = Vec::new();
        for hit in result.hits() {
            hits.push(PyHmmsearchHit {
                name: hit.name(),
                score: hit.score(),
                evalue: hit.evalue(),
            });
        }
        Ok(hits)
    }

    /// Run hmmsearch against a FASTA file on disk.
    fn search_fasta_file(&mut self, path: &str, hmm: &PyHmm) -> PyResult<Vec<PyHmmsearchHit>> {
        let result = self
            .inner
            .run_hmm_on_file(&hmm.inner, std::path::Path::new(path));
        let mut hits = Vec::new();
        for hit in result.hits() {
            hits.push(PyHmmsearchHit {
                name: hit.name(),
                score: hit.score(),
                evalue: hit.evalue(),
            });
        }
        Ok(hits)
    }
}

#[pyclass(name = "HmmerAlign", unsendable)]
pub struct PyHmmerAlign {
    inner: HmmerAlign,
}

#[pymethods]
impl PyHmmerAlign {
    #[new]
    fn new(hmm: &PyHmm) -> Self {
        PyHmmerAlign {
            inner: HmmerAlign::new(&hmm.inner),
        }
    }

    /// Align sequences and return Stockholm-format string.
    fn align_sequences(&self, sequences: Vec<(String, String)>) -> PyResult<String> {
        let easel_seqs = tuples_to_easel_sequences(sequences)?;
        self.inner
            .align_sequences_into_stockholm(&easel_seqs)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    }
}

#[pymodule]
fn hmmer_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyHmm>()?;
    m.add_class::<PyHmmerPipeline>()?;
    m.add_class::<PyHmmsearchHit>()?;
    m.add_class::<PyHmmerAlign>()?;
    Ok(())
}
