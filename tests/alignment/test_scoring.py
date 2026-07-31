import pytest
from Bio.Align import substitution_matrices

from bioseq.alignment.scoring import score_pair


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


def test_named_matrix_unsupported_residue_diagnostic_includes_context_and_alphabet():
    with pytest.raises(ValueError) as error:
        score_pair("a", "?", loaded_matrix="BLOSUM62")

    message = str(error.value)
    assert "('A', '?')" in message
    assert "BLOSUM62" in message
    assert "supported symbols" in message.lower()
    for symbol in ("A", "X", "*"):
        assert symbol in message


def test_raw_biopython_matrix_unsupported_residue_uses_generic_diagnostic():
    raw_matrix = substitution_matrices.load("BLOSUM62")

    with pytest.raises(ValueError) as error:
        score_pair("A", "?", loaded_matrix=raw_matrix)

    message = str(error.value)
    assert "('A', '?')" in message
    assert "selected substitution matrix" in message.lower()
    assert "supported symbols" in message.lower()
    assert "X" in message
    assert "*" in message
    assert "BLOSUM62" not in message
    assert "Array([" not in message


def test_raw_biopython_matrix_valid_scoring_is_unchanged():
    raw_matrix = substitution_matrices.load("BLOSUM62")

    assert score_pair("H", "P", loaded_matrix=raw_matrix) == -2.0
