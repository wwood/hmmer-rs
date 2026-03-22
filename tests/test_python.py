"""Tests for hmmer_rs Python bindings (PyO3)."""

import hmmer_rs


def test_read_hmm():
    hmms = hmmer_rs.Hmm.read_hmms_from_path(
        "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm"
    )
    assert len(hmms) == 1
    hmm = hmms[0]
    assert isinstance(hmm.name, str)
    assert len(hmm.name) > 0
    assert hmm.length > 0


def test_hmmsearch_sequences():
    hmms = hmmer_rs.Hmm.read_hmms_from_path(
        "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm"
    )
    hmm = hmms[0]

    pipeline = hmmer_rs.HmmerPipeline(hmm)
    sequences = [
        (
            "seq1",
            "MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
        ),
        ("seq2_no_hit", "AAAAAAAAAAAAAAAAAAAA"),
    ]
    hits = pipeline.search_sequences(sequences)

    assert len(hits) == 1
    hit = hits[0]
    assert hit.name == "seq1"
    assert hit.score > 100.0
    assert hit.evalue < 1e-40


def test_hmmsearch_with_evalue_threshold():
    hmms = hmmer_rs.Hmm.read_hmms_from_path(
        "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm"
    )
    hmm = hmms[0]

    # Very strict threshold — should still pass (evalue ~1.5e-48)
    pipeline = hmmer_rs.HmmerPipeline(hmm)
    pipeline.with_dom_evalue(1e-40)
    sequences = [
        (
            "good_hit",
            "MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
        ),
    ]
    hits = pipeline.search_sequences(sequences)
    assert len(hits) == 1

    # Impossibly strict — nothing should pass
    pipeline2 = hmmer_rs.HmmerPipeline(hmm)
    pipeline2.with_seq_evalue(1e-200)
    hits2 = pipeline2.search_sequences(sequences)
    assert len(hits2) == 0


def test_hmmalign():
    hmms = hmmer_rs.Hmm.read_hmms_from_path(
        "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm"
    )
    hmm = hmms[0]

    aligner = hmmer_rs.HmmerAlign(hmm)
    sequences = [
        (
            "test_seq",
            "MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
        ),
    ]
    stockholm = aligner.align_sequences(sequences)

    assert "# STOCKHOLM 1.0" in stockholm
    assert "test_seq" in stockholm
    assert "//" in stockholm


def test_hmmalign_multiple_sequences():
    hmms = hmmer_rs.Hmm.read_hmms_from_path(
        "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm"
    )
    hmm = hmms[0]

    aligner = hmmer_rs.HmmerAlign(hmm)
    sequences = [
        (
            "seq_one",
            "MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
        ),
        (
            "seq_two",
            "MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
        ),
    ]
    stockholm = aligner.align_sequences(sequences)

    assert "# STOCKHOLM 1.0" in stockholm
    assert "seq_one" in stockholm
    assert "seq_two" in stockholm
    assert "//" in stockholm


def test_pipeline_reuse():
    """Test that a pipeline can be used for multiple search calls."""
    hmms = hmmer_rs.Hmm.read_hmms_from_path(
        "tests/data/DNGNGWU00010_mingle_output_good_seqs.hmm"
    )
    hmm = hmms[0]

    pipeline = hmmer_rs.HmmerPipeline(hmm)
    sequences = [
        (
            "seq1",
            "MVYSGPNAPIEVGNSLPLSEIPLATEIHNIELTPGKGGQLVRSAGSSAQLLAKEGNYVTLRLPSGEMRFVRKECYATIGQ",
        ),
    ]

    hits1 = pipeline.search_sequences(sequences)
    assert len(hits1) == 1

    hits2 = pipeline.search_sequences(sequences)
    assert len(hits2) == 1


if __name__ == "__main__":
    test_read_hmm()
    test_hmmsearch_sequences()
    test_hmmsearch_with_evalue_threshold()
    test_hmmalign()
    test_hmmalign_multiple_sequences()
    test_pipeline_reuse()
    print("All Python tests passed!")
