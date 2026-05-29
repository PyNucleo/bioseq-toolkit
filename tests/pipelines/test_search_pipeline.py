from database.database_utils import normalize_database
from bioseq.pipelines.search_pipeline import search


def make_search_db():
    return normalize_database([
        "ATGCGT",
        "ATGCGA",
        "ATGAAA",
        "GGGGGG"
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
    db = normalize_database([
        "ATGAAA",
        "ATGCGT",
        "GGGGGG"
    ])

    results = search(
        query="ATGCG",
        database=db,
        k=3,
        threshold=1
    )

    assert results[0]["sequence"] == "ATGCGT"
    assert results[0]["shared_kmers"] == 3


def test_threshold_filtering():
    results = search(
        query="ATGCG",
        database=make_search_db(),
        k=3,
        threshold=2
    )

    assert [hit["sequence"] for hit in results] == ["ATGCGT", "ATGCGA"]
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
    assert all(hit["shared_kmers"] == 3 for hit in results)


def test_empty_results():
    db = normalize_database([
        "AAAAAA",
        "CCCCCC"
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