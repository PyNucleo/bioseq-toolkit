import json
import subprocess
import sys

import pytest
from hypothesis import given, settings, strategies as st

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


@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize(
    ("k", "threshold", "exception", "message"),
    [
        (0, 1, ValueError, "k.*positive"),
        (-1, 1, ValueError, "k.*positive"),
        (1.5, 1, TypeError, "k.*integer"),
        (True, 1, TypeError, "k.*integer"),
        (2, 0, ValueError, "threshold.*positive"),
        (2, -1, ValueError, "threshold.*positive"),
        (2, 1.5, TypeError, "threshold.*integer"),
        (2, True, TypeError, "threshold.*integer"),
    ],
)
def test_multi_search_validates_k_and_threshold_before_selecting_mode(
    indexed, k, threshold, exception, message
):
    with pytest.raises(exception, match=message):
        multi_search(
            [{"id": "query", "sequence": "ATGC"}],
            make_db([{"id": "seq", "sequence": "ATGC"}]),
            k=k,
            threshold=threshold,
            indexed=indexed,
        )


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

    Regular multi-search calls search(), whose k-mer candidate stage applies
    the fixed 0.3 ratio filter before rank_by_shared_kmers().

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


def test_indexed_multi_search_refinement_uses_real_simple_scoring_and_preserves_schema():
    db = make_db([
        {"id": "seq1", "sequence": "AAAA"},
    ])

    queries = [
        {"id": "query_seq", "sequence": "AAAA"},
    ]

    results = multi_search(
        queries,
        database=db,
        k=1,
        threshold=1,
        top_n_hits=1,
        indexed=True,
        refinement=True,
        match_score=3,
        mismatch_score=-7,
        gap_penalty=-2,
        matrix=None,
    )

    assert set(results[0]) == {"query_id", "query_sequence", "query_hits"}
    assert results[0]["query_id"] == "query_seq"
    assert results[0]["query_sequence"] == "AAAA"

    hit = results[0]["query_hits"][0]
    assert hit["id"] == "seq1"
    assert hit["sequence"] == "AAAA"
    assert hit["shared_kmers"] == 1
    assert hit["sw_score"] == 12
    assert hit["best_positions"] == [(4, 4)]


def test_indexed_multi_search_handles_empty_and_short_inputs():
    database = make_db([
        {"id": "short", "sequence": "AT"},
        {"id": "long", "sequence": "ATGCG"},
    ])

    assert multi_search([], database, k=3, threshold=1, indexed=True) == []

    results = multi_search(
        [
            {"id": "short_query", "sequence": "AT"},
            {"id": "normal_query", "sequence": "ATG"},
        ],
        database,
        k=3,
        threshold=1,
        indexed=True,
    )

    assert results[0]["query_hits"] == []
    assert simplify_hits(results[1]["query_hits"]) == [
        {"id": "long", "shared_kmers": 1}
    ]

    empty_database_results = multi_search(
        [{"id": "query", "sequence": "ATG"}],
        make_db([]),
        k=3,
        threshold=1,
        indexed=True,
    )

    assert empty_database_results == [{
        "query_id": "query",
        "query_sequence": "ATG",
        "query_hits": [],
    }]


def test_regular_and_indexed_multi_search_match_for_mixed_case_inputs():
    database = make_db([
        {"id": "seq_b", "sequence": "atgcgt"},
        {"id": "seq_a", "sequence": "ATGCGT"},
    ])
    queries = [{"id": "mixed", "sequence": "aTgCg"}]

    regular_results = multi_search(
        queries, database, k=3, threshold=1, indexed=False
    )
    indexed_results = multi_search(
        queries, database, k=3, threshold=1, indexed=True
    )

    assert indexed_results == regular_results
    assert indexed_results[0]["query_sequence"] == "aTgCg"
    assert simplify_hits(indexed_results[0]["query_hits"]) == [
        {"id": "seq_a", "shared_kmers": 3},
        {"id": "seq_b", "shared_kmers": 3},
    ]


def test_indexed_multi_search_rejects_duplicate_database_ids():
    database = make_db([
        {"id": "duplicate", "sequence": "ATGC"},
        {"id": "unique", "sequence": "GGGG"},
        {"id": "duplicate", "sequence": "CCCC"},
    ])

    with pytest.raises(ValueError) as error:
        multi_search(
            [{"id": "query", "sequence": "ATGC"}],
            database,
            k=2,
            threshold=1,
            indexed=True,
        )

    message = str(error.value)
    assert "Duplicate database IDs detected:" in message
    assert "'duplicate' appears in records [1, 3]" in message


def test_regular_and_indexed_multi_search_are_equivalent_for_multiple_queries():
    database = make_db([
        {"id": "z_tie", "sequence": "MKTAAA"},
        {"id": "a_tie", "sequence": "MKTAAA"},
        {"id": "weak_alpha", "sequence": "MKTGGG"},
        {"id": "beta", "sequence": "GGGCCC"},
        {"id": "below_threshold", "sequence": "GGGAAA"},
    ])
    queries = [
        {"id": "alpha", "sequence": "MKTAAA"},
        {"id": "beta", "sequence": "GGGCCC"},
    ]

    regular_results = multi_search(
        queries, database, k=3, threshold=2, top_n_hits=2, indexed=False
    )
    indexed_results = multi_search(
        queries, database, k=3, threshold=2, top_n_hits=2, indexed=True
    )

    assert indexed_results == regular_results
    assert simplify_hits(indexed_results[0]["query_hits"]) == [
        {"id": "a_tie", "shared_kmers": 4},
        {"id": "z_tie", "shared_kmers": 4},
    ]
    assert simplify_hits(indexed_results[1]["query_hits"]) == [
        {"id": "beta", "shared_kmers": 4}
    ]


def _query_record_strategy():
    return st.fixed_dictionaries(
        {
            "id": st.text(min_size=3, max_size=10),
            "sequence": st.text(
                alphabet="ATGCatgc", min_size=4, max_size=10
            ),
        }
    )


@settings(max_examples=150, deadline=None)
@given(
    query=st.lists(
        _query_record_strategy(),
        min_size=1,
        max_size=10,
        unique_by=lambda record: record["id"],
    ),
    database=st.lists(
        st.text(alphabet="ATGCatgc", min_size=4, max_size=10),
        min_size=1,
        max_size=10,
    ),
    k=st.integers(min_value=1, max_value=8),
    threshold=st.integers(min_value=1, max_value=4),
    top_n_hits=st.integers(min_value=1, max_value=4),
)
def test_regular_and_indexed_multi_search_are_equivalent_for_multiple_queries_with_hypothesis(
    query, database, k, threshold, top_n_hits
):
    regular_multi_search = multi_search(
        query,
        database,
        k,
        threshold,
        top_n_hits,
        indexed=False,
        refinement=False,
    )
    indexed_multi_search = multi_search(
        query,
        database,
        k,
        threshold,
        top_n_hits,
        indexed=True,
        refinement=False,
    )

    assert regular_multi_search == indexed_multi_search
