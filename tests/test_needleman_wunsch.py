from bioseq.alignment.needleman_wunsch import global_alignment


def test_global_alignment_identical():
    result = global_alignment("GATTACA", "GATTACA")

    expected = [("GATTACA", "GATTACA")]

    assert result == expected


def test_global_alignment_single_gap():
    result = global_alignment("ACTG", "ACG")

    assert len(result) > 0