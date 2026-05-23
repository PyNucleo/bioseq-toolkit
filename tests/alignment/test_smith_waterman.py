from bioseq.alignment.smith_waterman import local_alignment


def test_local_alignment_basic():
    result = local_alignment("GATTACA", "TTAC")

    assert len(result) > 0


def test_local_alignment_identical():
    result = local_alignment("AAAA", "AAAA")

    assert len(result) > 0