from database.sequence_database import SequenceDatabase
from bioseq.search.kmer_search import filter_by_relative_score, kmer_search


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