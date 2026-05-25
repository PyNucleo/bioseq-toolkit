# tests/search/test_refinement.py

from database.sequence_database import SequenceDatabase

from bioseq.search.refinement import refine_hits
from bioseq.pipelines.search_pipeline import search


def test_refinement_adds_sw_score_and_best_positions():
    top_hits = [
        {
            "seq_id": None,
            "sequence": "ATGCGT",
            "shared_kmers": 3
        }
    ]

    results = refine_hits(query="ATGCG", top_hits=top_hits)

    assert len(results) == 1
    assert results[0]["sequence"] == "ATGCGT"
    assert results[0]["sw_score"] == 10
    assert results[0]["best_positions"] == [(5, 5)]


def test_refinement_ranks_by_sw_score_not_shared_kmers():
    top_hits = [
        {
            "seq_id": None,
            "sequence": "ATGAAA",
            "shared_kmers": 99
        },
        {
            "seq_id": None,
            "sequence": "ATGCGT",
            "shared_kmers": 1
        }
    ]

    results = refine_hits(query="ATGCG", top_hits=top_hits)

    assert results[0]["sequence"] == "ATGCGT"
    assert results[0]["sw_score"] >= results[1]["sw_score"]


def test_refinement_handles_empty_hits():
    results = refine_hits(query="ATGCG", top_hits=[])

    assert results == []


def test_search_pipeline_with_refinement_adds_sw_scores():
    db = SequenceDatabase([
        "ATGAAA",
        "ATGCGT",
        "GGGGGG"
    ])

    results = search(
        query="ATGCG",
        database=db,
        k=3,
        threshold=1,
        refinement=True
    )

    assert len(results) > 0
    assert "sw_score" in results[0]
    assert "best_positions" in results[0]
    assert results[0]["sequence"] == "ATGCGT"