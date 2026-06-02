from bioseq.alignment.scoring import score_pair


def test_blosum62_known_pair_scores():
    assert score_pair("H", "H", "BLOSUM62") == 8.0
    assert score_pair("E", "E", "BLOSUM62") == 5.0
    assert score_pair("H", "P", "BLOSUM62") == -2.0
