from bioseq.search.kmer_search import (
    generate_kmers,
    get_shared_kmers,
    filter_candidates_ratio,
    kmer_search,
)
from database.sequence_database import SequenceDatabase

def test_generate_kmers():

    result = generate_kmers("GATTACA", 3)

    expected = {
        "GAT",
        "ATT",
        "TTA",
        "TAC",
        "ACA",
    }

    assert result == expected


def test_get_shared_kmers_identical():

    temp_query = generate_kmers("GATTACA", 3)

    result = get_shared_kmers(
        temp_query,
        "GATTACA",
        3
    )

    assert result == 5


def test_get_shared_kmers_none():

    result = get_shared_kmers(
        "AAAAAAA",
        "CCCCCCC",
        3
    )

    assert result == 0


def test_filter_candidates_ratio():

    candidates = {
        "seq1": 10,
        "seq2": 8,
        "seq3": 1,
    }

    result = filter_candidates_ratio(
        10,
        candidates,
        0.5
    )

    expected = {
        "seq1": 10,
        "seq2": 8,
    }

    assert result == expected


def test_kmer_search_basic():

    database = SequenceDatabase([
        "GATTACA",
        "CCCCCCC",
        "GATTTTT",
    ])

    result = kmer_search(
        "GATTACA",
        database,
        3,
        1
    )

    assert "GATTACA" in result
    assert "GATTTTT" in result


def test_generate_kmers_large_k():

    result = generate_kmers("ATG", 5)

    expected = set()

    assert result == expected