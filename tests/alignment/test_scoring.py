import pytest

from bioseq.alignment.scoring import score_pair
from bioseq.alignment.substitution_matrices import load_matrix


def test_blosum62_known_pair_scores():
    assert score_pair("H", "H", "BLOSUM62") == 8.0
    assert score_pair("E", "E", "BLOSUM62") == 5.0
    assert score_pair("H", "P", "BLOSUM62") == -2.0


def test_blosum62_normalizes_lowercase_residues():
    assert score_pair("h", "p", loaded_matrix="BLOSUM62") == -2.0


def test_simple_scoring_allows_unrestricted_characters():
    assert score_pair("A", "A", loaded_matrix=None, match_score=3, mismatch_score=-4) == 3
    assert score_pair("A", "?", loaded_matrix=None, match_score=3, mismatch_score=-4) == -4


@pytest.mark.parametrize(
    ("a", "b", "normalized_a", "normalized_b"),
    [("?", "A", "?", "A"), ("A", "?", "A", "?")],
)
def test_blosum62_rejects_unsupported_residues(a, b, normalized_a, normalized_b):
    with pytest.raises(ValueError) as error:
        score_pair(a, b, loaded_matrix="BLOSUM62")

    message = str(error.value).lower()
    assert "not supported" in message
    assert "substitution matrix" in message
    assert repr(normalized_a) in str(error.value)
    assert repr(normalized_b) in str(error.value)


def test_loaded_biopython_matrix_is_supported_and_validated():
    matrix = load_matrix("BLOSUM62")

    assert score_pair("H", "P", loaded_matrix=matrix) == -2.0

    with pytest.raises(ValueError, match="not supported"):
        score_pair("H", "?", loaded_matrix=matrix)
