import pytest

from bioseq.alignment import substitution_matrices as project_matrices
from bioseq.alignment.substitution_matrices import load_matrix


def test_load_matrix_translates_unknown_name_error():
    load_matrix.cache_clear()

    with pytest.raises(ValueError) as error:
        load_matrix("BLOSSUM62")

    assert "unknown substitution matrix" in str(error.value).lower()
    assert "BLOSSUM62" in str(error.value)
    assert isinstance(error.value.__cause__, FileNotFoundError)


def test_load_matrix_does_not_translate_unexpected_loader_errors(monkeypatch):
    def raise_unexpected_error(matrix_name):
        raise RuntimeError("loader failed")

    load_matrix.cache_clear()
    monkeypatch.setattr(project_matrices.substitution_matrices, "load", raise_unexpected_error)

    with pytest.raises(RuntimeError, match="loader failed"):
        load_matrix("BLOSUM62")

    load_matrix.cache_clear()
