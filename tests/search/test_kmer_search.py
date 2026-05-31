import pytest

from database.sequence_database import SequenceDatabase
from bioseq.search.kmer_search import (
    filter_by_relative_score,
    generate_kmers,
    kmer_search,
    validate_kmer_params,
)


def test_filter_by_relative_score():
    candidates = [
        {"id": "seq1", "sequence": "ATGCGT", "shared_kmers": 10},
        {"id": "seq2", "sequence": "ATGAAA", "shared_kmers": 3},
        {"id": "seq3", "sequence": "GGGGGG", "shared_kmers": 2},
        {"id": "seq4", "sequence": "CCCCCC", "shared_kmers": 0},
    ]

    filtered = filter_by_relative_score(
        max_kmers=10,
        candidates_dict=candidates,
        ratio=0.3
    )

    assert isinstance(filtered, list)
    assert [hit["id"] for hit in filtered] == ["seq1", "seq2"]
    assert filtered[0]["sequence"] == "ATGCGT"
    assert filtered[0]["shared_kmers"] == 10

    assert filter_by_relative_score(
        max_kmers=0,
        candidates_dict=[],
        ratio=0.3
    ) == []


def test_kmer_search_basic():
    db = SequenceDatabase([
        {"id": "seq1", "sequence": "ATGCGT"},
        {"id": "seq2", "sequence": "ATGCGT"},
        {"id": "seq3", "sequence": "ATGAAA"},
        {"id": "seq4", "sequence": "GGGGGG"},
    ])

    results = kmer_search(
        query="ATGCG",
        db=db,
        k=3,
        threshold=1
    )

    assert isinstance(results, list)
    assert all(isinstance(hit, dict) for hit in results)
    assert {"id", "sequence", "shared_kmers"}.issubset(results[0].keys())

    by_id = {hit["id"]: hit for hit in results}

    assert by_id["seq1"]["shared_kmers"] == 3
    assert by_id["seq2"]["shared_kmers"] == 3
    assert by_id["seq3"]["shared_kmers"] == 1

    assert "seq4" not in by_id

    assert by_id["seq1"]["sequence"] == by_id["seq2"]["sequence"]
    assert by_id["seq1"]["id"] != by_id["seq2"]["id"]

def test_validate_kmer_params_accepts_valid_inputs():
    assert validate_kmer_params("ATGC", k=2) is True
    assert validate_kmer_params("ATGC", k=2, threshold=1) is True


def test_validate_kmer_params_rejects_non_integer_k():
    with pytest.raises(TypeError):
        validate_kmer_params("ATGC", k=2.5)

    with pytest.raises(TypeError):
        validate_kmer_params("ATGC", k="2")


def test_validate_kmer_params_rejects_non_positive_k():
    with pytest.raises(ValueError):
        validate_kmer_params("ATGC", k=0)

    with pytest.raises(ValueError):
        validate_kmer_params("ATGC", k=-1)


def test_validate_kmer_params_returns_false_when_k_larger_than_sequence():
    assert validate_kmer_params("ATGC", k=10) is False


def test_validate_kmer_params_rejects_non_integer_threshold():
    with pytest.raises(TypeError):
        validate_kmer_params("ATGC", k=2, threshold=1.5)

    with pytest.raises(TypeError):
        validate_kmer_params("ATGC", k=2, threshold="1")


def test_validate_kmer_params_rejects_negative_threshold():
    with pytest.raises(ValueError):
        validate_kmer_params("ATGC", k=2, threshold=-1)


def test_generate_kmers_rejects_invalid_k_values():
    with pytest.raises(ValueError):
        generate_kmers("ATGC", 0)

    with pytest.raises(ValueError):
        generate_kmers("ATGC", -1)

    with pytest.raises(TypeError):
        generate_kmers("ATGC", 2.5)


def test_generate_kmers_returns_empty_set_when_k_larger_than_sequence():
    assert generate_kmers("ATGC", 10) == set()


def test_kmer_search_rejects_invalid_k_values():
    db = SequenceDatabase([
        {"id": "seq1", "sequence": "ATGCGT"},
    ])

    with pytest.raises(ValueError):
        kmer_search("ATGC", db, k=0, threshold=1)

    with pytest.raises(ValueError):
        kmer_search("ATGC", db, k=-1, threshold=1)

    with pytest.raises(TypeError):
        kmer_search("ATGC", db, k=2.5, threshold=1)


def test_kmer_search_rejects_invalid_threshold_values():
    db = SequenceDatabase([
        {"id": "seq1", "sequence": "ATGCGT"},
    ])

    with pytest.raises(ValueError):
        kmer_search("ATGC", db, k=3, threshold=-1)

    with pytest.raises(TypeError):
        kmer_search("ATGC", db, k=3, threshold=1.5)


def test_kmer_search_returns_empty_when_k_larger_than_query():
    db = SequenceDatabase([
        {"id": "seq1", "sequence": "ATGCGT"},
    ])

    assert kmer_search("ATGC", db, k=10, threshold=1) == []