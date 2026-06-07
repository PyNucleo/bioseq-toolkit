import json
import subprocess
import sys

import pytest

from database.sequence_database import SequenceDatabase
from bioseq.pipelines.search_pipeline import multi_search


def make_db(records):
    return SequenceDatabase(records)


def simplify_hits(hits):
    return [
        {
            "id": hit["id"],
            "shared_kmers": hit["shared_kmers"],
        }
        for hit in hits
    ]


def test_indexed_multi_search_sorts_before_top_n():
    """
    Regression test for the exact bug you found.

    The best hit is intentionally placed late in the database.
    If top_n is applied before sorting, this fails.
    """

    db = make_db([
        {"id": "weak_1", "sequence": "MKTGGG"},
        {"id": "weak_2", "sequence": "MKTPPP"},
        {"id": "best_late", "sequence": "MKTAAA"},
    ])

    queries = [
        {"id": "query_seq", "sequence": "MKTAAA"},
    ]

    results = multi_search(
        queries,
        database=db,
        k=3,
        threshold=1,
        top_n_hits=1,
        indexed=True,
        refinement=False,
    )

    hits = results[0]["query_hits"]

    assert simplify_hits(hits) == [
        {"id": "best_late", "shared_kmers": 4}
    ]


def test_indexed_multi_search_uses_deterministic_id_tie_breaking():
    """
    Equal shared_kmers means biologically tied under this scoring rule.
    ID tie-breaking only makes output deterministic.
    """

    db = make_db([
        {"id": "seq_c", "sequence": "MKTAAA"},
        {"id": "seq_a", "sequence": "MKTAAA"},
        {"id": "seq_b", "sequence": "MKTAAA"},
    ])

    queries = [
        {"id": "query_seq", "sequence": "MKTAAA"},
    ]

    results = multi_search(
        queries,
        database=db,
        k=3,
        threshold=1,
        top_n_hits=3,
        indexed=True,
        refinement=False,
    )

    hits = results[0]["query_hits"]

    assert simplify_hits(hits) == [
        {"id": "seq_a", "shared_kmers": 4},
        {"id": "seq_b", "shared_kmers": 4},
        {"id": "seq_c", "shared_kmers": 4},
    ]


def test_indexed_multi_search_handles_multiple_queries_independently():
    db = make_db([
        {"id": "alpha_exact", "sequence": "MKTAAA"},
        {"id": "alpha_partial", "sequence": "MKTAAG"},
        {"id": "beta_exact", "sequence": "GGGCCC"},
        {"id": "beta_partial", "sequence": "GGGCCA"},
        {"id": "unrelated", "sequence": "PPPPPP"},
    ])

    queries = [
        {"id": "query_alpha", "sequence": "MKTAAA"},
        {"id": "query_beta", "sequence": "GGGCCC"},
    ]

    results = multi_search(
        queries,
        database=db,
        k=3,
        threshold=1,
        top_n_hits=2,
        indexed=True,
        refinement=False,
    )

    assert [result["query_id"] for result in results] == [
        "query_alpha",
        "query_beta",
    ]

    assert simplify_hits(results[0]["query_hits"]) == [
        {"id": "alpha_exact", "shared_kmers": 4},
        {"id": "alpha_partial", "shared_kmers": 3},
    ]

    assert simplify_hits(results[1]["query_hits"]) == [
        {"id": "beta_exact", "shared_kmers": 4},
        {"id": "beta_partial", "shared_kmers": 3},
    ]


def test_regular_multi_search_uses_query_hits_key():
    db = make_db([
        {"id": "seq1", "sequence": "MKTAAA"},
        {"id": "seq2", "sequence": "MKTGGG"},
    ])

    queries = [
        {"id": "query_seq", "sequence": "MKTAAA"},
    ]

    results = multi_search(
        queries,
        database=db,
        k=3,
        threshold=1,
        top_n_hits=10,
        indexed=False,
        refinement=False,
    )

    assert "query_hits" in results[0]
    assert "indexed_hits" not in results[0]

def test_regular_and_indexed_multi_search_use_same_output_key():
    db = make_db([
        {"id": "seq1", "sequence": "MKTAAA"},
    ])

    queries = [
        {"id": "query_seq", "sequence": "MKTAAA"},
    ]

    regular_results = multi_search(
        queries,
        database=db,
        k=3,
        threshold=1,
        top_n_hits=1,
        indexed=False,
        refinement=False,
    )

    indexed_results = multi_search(
        queries,
        database=db,
        k=3,
        threshold=1,
        top_n_hits=1,
        indexed=True,
        refinement=False,
    )

    assert "query_hits" in regular_results[0]
    assert "query_hits" in indexed_results[0]
    assert "indexed_hits" not in regular_results[0]
    assert "indexed_hits" not in indexed_results[0]

def test_regular_multi_search_applies_search_pipeline_ratio_filter():
    """
    This documents the current behavior.

    Regular multi-search calls search(), and search() uses rank_by_shared_kmers().
    If rank_by_shared_kmers applies the 0.3 ratio filter, weak seq2 is removed.

    Query k-mers:
    MKT, KTA, TAA, AAA

    seq1 shares 4.
    seq2 shares 1.

    Since 1 < 0.3 * 4, seq2 should be filtered out.
    """

    db = make_db([
        {"id": "seq1", "sequence": "MKTAAA"},
        {"id": "seq2", "sequence": "MKTGGG"},
    ])

    queries = [
        {"id": "query_seq", "sequence": "MKTAAA"},
    ]

    results = multi_search(
        queries,
        database=db,
        k=3,
        threshold=1,
        top_n_hits=10,
        indexed=False,
        refinement=False,
    )

    assert simplify_hits(results[0]["query_hits"]) == [
        {"id": "seq1", "shared_kmers": 4}
    ]


def test_regular_and_indexed_multi_search_apply_same_ratio_filter():
    """
    Regular and indexed multi-search should now be semantically equivalent.

    Query k-mers for MKTAAA with k=3:
    MKT, KTA, TAA, AAA

    seq1 shares 4.
    seq2 shares 1.

    With the 0.3 ratio filter:
    max_shared_kmers = 4
    cutoff = 0.3 * 4 = 1.2

    Therefore seq2 should be removed in both modes.
    """

    db = make_db([
        {"id": "seq1", "sequence": "MKTAAA"},
        {"id": "seq2", "sequence": "MKTGGG"},
    ])

    queries = [
        {"id": "query_seq", "sequence": "MKTAAA"},
    ]

    regular_results = multi_search(
        queries,
        database=db,
        k=3,
        threshold=1,
        top_n_hits=10,
        indexed=False,
        refinement=False,
    )

    indexed_results = multi_search(
        queries,
        database=db,
        k=3,
        threshold=1,
        top_n_hits=10,
        indexed=True,
        refinement=False,
    )

    assert simplify_hits(regular_results[0]["query_hits"]) == [
        {"id": "seq1", "shared_kmers": 4}
    ]

    assert simplify_hits(indexed_results[0]["query_hits"]) == [
        {"id": "seq1", "shared_kmers": 4}
    ]

    assert simplify_hits(indexed_results[0]["query_hits"]) == simplify_hits(
        regular_results[0]["query_hits"]
    )

def test_multi_search_accepts_query_records_directly():
    db = make_db([
        {"id": "seq1", "sequence": "MKTAAA"},
    ])

    queries = [
        {"id": "query_seq", "sequence": "MKTAAA"},
    ]

    results = multi_search(
        queries,
        database=db,
        k=3,
        threshold=1,
        top_n_hits=1,
        indexed=False,
        refinement=False,
    )

    assert results[0]["query_id"] == "query_seq"
    assert results[0]["query_sequence"] == "MKTAAA"
    assert results[0]["query_hits"][0]["id"] == "seq1"


def test_multi_search_accepts_query_fasta_path(tmp_path):
    query_file = tmp_path / "queries.fasta"

    query_file.write_text(
        ">query_seq\n"
        "MKTAAA\n"
    )

    db = make_db([
        {"id": "seq1", "sequence": "MKTAAA"},
    ])

    results = multi_search(
        query_file,
        database=db,
        k=3,
        threshold=1,
        top_n_hits=1,
        indexed=False,
        refinement=False,
    )

    assert results[0]["query_id"] == "query_seq"
    assert results[0]["query_sequence"] == "MKTAAA"
    assert results[0]["query_hits"][0]["id"] == "seq1"


def test_indexed_multi_search_refinement_is_not_silently_ignored():
    """
    refinement=True with indexed=True should not be silently ignored.
    Either implement indexed refinement or raise NotImplementedError.
    """

    db = make_db([
        {"id": "seq1", "sequence": "MKTAAA"},
    ])

    queries = [
        {"id": "query_seq", "sequence": "MKTAAA"},
    ]

    with pytest.raises(NotImplementedError):
        multi_search(
            queries,
            database=db,
            k=3,
            threshold=1,
            top_n_hits=1,
            indexed=True,
            refinement=True,
        )