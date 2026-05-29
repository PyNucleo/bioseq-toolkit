from database.sequence_database import SequenceDatabase
from bioseq.pipelines.search_pipeline import search


def make_search_db():
    return SequenceDatabase([
        {"id": "seq1", "sequence": "ATGCGT"},
        {"id": "seq2", "sequence": "ATGCGA"},
        {"id": "seq3", "sequence": "ATGAAA"},
        {"id": "seq4", "sequence": "GGGGGG"},
    ])


def test_pipeline_returns_hits():
    results = search(
        query="ATGCG",
        database=make_search_db(),
        k=3,
        threshold=1
    )

    assert len(results) > 0
    assert all(isinstance(hit, dict) for hit in results)
    assert {"id", "sequence", "shared_kmers"}.issubset(results[0].keys())


def test_best_hit_ranked_first():
    db = SequenceDatabase([
        {"id": "weak", "sequence": "ATGAAA"},
        {"id": "best", "sequence": "ATGCGT"},
        {"id": "none", "sequence": "GGGGGG"},
    ])

    results = search(
        query="ATGCG",
        database=db,
        k=3,
        threshold=1
    )

    assert results[0]["id"] == "best"
    assert results[0]["sequence"] == "ATGCGT"
    assert results[0]["shared_kmers"] == 3


def test_threshold_filtering():
    results = search(
        query="ATGCG",
        database=make_search_db(),
        k=3,
        threshold=2
    )

    assert [hit["id"] for hit in results] == ["seq1", "seq2"]
    assert all(hit["shared_kmers"] >= 2 for hit in results)


def test_top_n_hits_filtering():
    results = search(
        query="ATGCG",
        database=make_search_db(),
        k=3,
        threshold=1,
        top_n_hits=2
    )

    assert len(results) == 2
    assert [hit["id"] for hit in results] == ["seq1", "seq2"]
    assert all(hit["shared_kmers"] == 3 for hit in results)


def test_empty_results():
    db = SequenceDatabase([
        {"id": "seq1", "sequence": "AAAAAA"},
        {"id": "seq2", "sequence": "CCCCCC"},
    ])

    results = search(
        query="ATGCG",
        database=db,
        k=3,
        threshold=1
    )

    assert results == []


def test_search_pipeline_with_refinement_adds_sw_scores():
    results = search(
        query="ATGCG",
        database=make_search_db(),
        k=3,
        threshold=1,
        top_n_hits=3,
        refinement=True
    )

    assert len(results) == 3

    required_keys = {
        "id",
        "sequence",
        "shared_kmers",
        "sw_score",
        "best_positions"
    }

    assert all(required_keys.issubset(hit.keys()) for hit in results)
    assert all(isinstance(hit["sw_score"], int) for hit in results)
    assert all(isinstance(hit["best_positions"], list) for hit in results)

    assert results == sorted(
        results,
        key=lambda hit: hit["sw_score"],
        reverse=True
    )

    assert results[0]["sw_score"] >= results[-1]["sw_score"]