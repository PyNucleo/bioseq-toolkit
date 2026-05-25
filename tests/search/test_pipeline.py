# tests/pipelines/test_search_pipeline.py

from database.sequence_database import SequenceDatabase
from bioseq.pipelines.search_pipeline import search


def test_pipeline_returns_hits():
    db = SequenceDatabase([
        "ATGCGT",
        "ATGAAA",
        "GGGGGG",
        "ATGCGA"
    ])

    results = search(query="ATGCG", database=db)

    assert len(results) > 0


def test_best_hit_ranked_first():
    db = SequenceDatabase([
        "ATGAAA",
        "ATGCGT",
        "GGGGGG"
    ])

    results = search(query="ATGCG", database=db, k=3)

    assert results[0]["sequence"] == "ATGCGT"


def test_threshold_filtering():
    db = SequenceDatabase([
        "ATGCGT",
        "ATGAAA",
        "GGGGGG"
    ])

    results = search(
        query="ATGCG",
        database=db,
        k=3,
        threshold=2
    )

    for hit in results:
        assert hit["shared_kmers"] >= 2


def test_top_n_hits_filtering():
    db = SequenceDatabase([
        "ATGCGT",
        "ATGCGA",
        "ATGAAA",
        "GGGGGG"
    ])

    results = search(
        query="ATGCG",
        database=db,
        top_n_hits=2
    )

    assert len(results) == 2


def test_empty_results():
    db = SequenceDatabase([
        "AAAAAA",
        "CCCCCC",
        "GGGGGG"
    ])

    results = search(
        query="TTTTTT",
        database=db,
        threshold=1
    )

    assert results == []