import pytest

from bioseq.pipelines.search_pipeline import multi_search
from database.sequence_database import SequenceDatabase


def make_tiny_database():
    return SequenceDatabase([
        {"id": "seq1", "sequence": "ATGCGT"},
        {"id": "seq2", "sequence": "ATGCGA"},
        {"id": "seq3", "sequence": "GGGGGG"},
        {"id": "seq4", "sequence": "TTTAAA"},
        {"id": "seq5", "sequence": "MKWVTFISLL"},
    ])


def write_query_fasta(tmp_path):
    query_fasta = tmp_path / "queries.fasta"

    query_fasta.write_text(
        ">query_exact_seq1_seq2\n"
        "ATGCG\n"
        ">query_g_rich\n"
        "GGGG\n"
        ">query_ta_region\n"
        "TTTAA\n"
        ">query_protein_like\n"
        "MKWV\n"
        ">query_no_hit\n"
        "CCCC\n"
    )

    return query_fasta


def hits_to_counts(hits):
    return {
        hit["id"]: hit["shared_kmers"]
        for hit in hits
    }


def simplify_hits(hits):
    return [
        {
            "id": hit["id"],
            "shared_kmers": hit["shared_kmers"],
        }
        for hit in hits
    ]


def result_by_query_id(results):
    return {
        result["query_id"]: result
        for result in results
    }


def test_multi_search_indexed_returns_one_result_per_query(tmp_path):
    database = make_tiny_database()
    query_fasta = write_query_fasta(tmp_path)

    results = multi_search(
        query_fasta=query_fasta,
        database=database,
        k=3,
        threshold=1,
        indexed=True,
        refinement=False,
    )

    assert len(results) == 5

    query_ids = [result["query_id"] for result in results]

    assert query_ids == [
        "query_exact_seq1_seq2",
        "query_g_rich",
        "query_ta_region",
        "query_protein_like",
        "query_no_hit",
    ]


def test_multi_search_indexed_uses_query_hits_key(tmp_path):
    database = make_tiny_database()
    query_fasta = write_query_fasta(tmp_path)

    results = multi_search(
        query_fasta=query_fasta,
        database=database,
        k=3,
        threshold=1,
        indexed=True,
        refinement=False,
    )

    for result in results:
        assert "query_hits" in result
        assert "indexed_hits" not in result


def test_multi_search_indexed_returns_expected_hits_and_counts(tmp_path):
    database = make_tiny_database()
    query_fasta = write_query_fasta(tmp_path)

    results = multi_search(
        query_fasta=query_fasta,
        database=database,
        k=3,
        threshold=1,
        indexed=True,
        refinement=False,
    )

    results_by_query = result_by_query_id(results)

    assert hits_to_counts(
        results_by_query["query_exact_seq1_seq2"]["query_hits"]
    ) == {
        "seq1": 3,
        "seq2": 3,
    }

    assert hits_to_counts(
        results_by_query["query_g_rich"]["query_hits"]
    ) == {
        "seq3": 1,
    }

    assert hits_to_counts(
        results_by_query["query_ta_region"]["query_hits"]
    ) == {
        "seq4": 3,
    }

    assert hits_to_counts(
        results_by_query["query_protein_like"]["query_hits"]
    ) == {
        "seq5": 2,
    }

    assert results_by_query["query_no_hit"]["query_hits"] == []


def test_multi_search_indexed_sorts_hits_by_shared_kmers_then_id(tmp_path):
    database = make_tiny_database()
    query_fasta = write_query_fasta(tmp_path)

    results = multi_search(
        query_fasta=query_fasta,
        database=database,
        k=3,
        threshold=1,
        indexed=True,
        refinement=False,
    )

    results_by_query = result_by_query_id(results)

    hits = results_by_query["query_exact_seq1_seq2"]["query_hits"]

    assert simplify_hits(hits) == [
        {"id": "seq1", "shared_kmers": 3},
        {"id": "seq2", "shared_kmers": 3},
    ]


def test_multi_search_indexed_applies_threshold_per_query(tmp_path):
    database = make_tiny_database()
    query_fasta = write_query_fasta(tmp_path)

    results = multi_search(
        query_fasta=query_fasta,
        database=database,
        k=3,
        threshold=3,
        indexed=True,
        refinement=False,
    )

    results_by_query = result_by_query_id(results)

    assert hits_to_counts(
        results_by_query["query_exact_seq1_seq2"]["query_hits"]
    ) == {
        "seq1": 3,
        "seq2": 3,
    }

    assert results_by_query["query_g_rich"]["query_hits"] == []

    assert hits_to_counts(
        results_by_query["query_ta_region"]["query_hits"]
    ) == {
        "seq4": 3,
    }

    assert results_by_query["query_protein_like"]["query_hits"] == []

    assert results_by_query["query_no_hit"]["query_hits"] == []


def test_multi_search_indexed_uses_database_wide_index_not_first_query_hits(tmp_path):
    database = make_tiny_database()

    query_fasta = tmp_path / "unrelated_queries.fasta"
    query_fasta.write_text(
        ">dna_query\n"
        "ATGCG\n"
        ">protein_query\n"
        "MKWV\n"
    )

    results = multi_search(
        query_fasta=query_fasta,
        database=database,
        k=3,
        threshold=1,
        indexed=True,
        refinement=False,
    )

    results_by_query = result_by_query_id(results)

    assert hits_to_counts(
        results_by_query["dna_query"]["query_hits"]
    ) == {
        "seq1": 3,
        "seq2": 3,
    }

    assert hits_to_counts(
        results_by_query["protein_query"]["query_hits"]
    ) == {
        "seq5": 2,
    }


def test_multi_search_indexed_respects_top_n_hits_per_query(tmp_path):
    database = make_tiny_database()
    query_fasta = write_query_fasta(tmp_path)

    results = multi_search(
        query_fasta=query_fasta,
        database=database,
        k=3,
        threshold=1,
        top_n_hits=1,
        indexed=True,
        refinement=False,
    )

    results_by_query = result_by_query_id(results)

    assert simplify_hits(
        results_by_query["query_exact_seq1_seq2"]["query_hits"]
    ) == [
        {"id": "seq1", "shared_kmers": 3}
    ]

    assert simplify_hits(
        results_by_query["query_g_rich"]["query_hits"]
    ) == [
        {"id": "seq3", "shared_kmers": 1}
    ]

    assert results_by_query["query_no_hit"]["query_hits"] == []


def test_regular_multi_search_propagates_simple_refinement_scores():
    results = multi_search(
        [{"id": "query", "sequence": "AAAA"}],
        SequenceDatabase([{"id": "exact", "sequence": "AAAA"}]),
        k=1,
        threshold=1,
        indexed=False,
        refinement=True,
        match_score=3,
        mismatch_score=-7,
        gap_penalty=-2,
        matrix=None,
    )

    assert results[0]["query_hits"][0]["sw_score"] == 12


def test_indexed_refinement_keeps_queries_independent_and_reranks_selected_hits():
    database = SequenceDatabase([
        {"id": "a_alpha_weak", "sequence": "AATA"},
        {"id": "z_alpha_exact", "sequence": "AAAA"},
        {"id": "a_beta_weak", "sequence": "CGCC"},
        {"id": "z_beta_exact", "sequence": "CCCC"},
    ])
    queries = [
        {"id": "alpha", "sequence": "AAAA"},
        {"id": "beta", "sequence": "CCCC"},
    ]

    unrefined = multi_search(
        queries,
        database,
        k=1,
        threshold=1,
        top_n_hits=2,
        indexed=True,
        refinement=False,
    )
    refined = multi_search(
        queries,
        database,
        k=1,
        threshold=1,
        top_n_hits=2,
        indexed=True,
        refinement=True,
        match_score=2,
        mismatch_score=-1,
        gap_penalty=-2,
    )

    assert [result["query_id"] for result in refined] == ["alpha", "beta"]
    assert [result["query_sequence"] for result in refined] == ["AAAA", "CCCC"]

    assert [hit["id"] for hit in unrefined[0]["query_hits"]] == [
        "a_alpha_weak",
        "z_alpha_exact",
    ]
    assert [hit["id"] for hit in refined[0]["query_hits"]] == [
        "z_alpha_exact",
        "a_alpha_weak",
    ]
    assert [hit["sw_score"] for hit in refined[0]["query_hits"]] == [8, 5]

    assert [hit["id"] for hit in unrefined[1]["query_hits"]] == [
        "a_beta_weak",
        "z_beta_exact",
    ]
    assert [hit["id"] for hit in refined[1]["query_hits"]] == [
        "z_beta_exact",
        "a_beta_weak",
    ]
    assert [hit["sw_score"] for hit in refined[1]["query_hits"]] == [8, 5]

    for result in refined:
        assert result["query_hits"] == sorted(
            result["query_hits"],
            key=lambda hit: hit["sw_score"],
            reverse=True,
        )
        assert {hit["id"] for hit in result["query_hits"]} == {
            hit["id"] for hit in unrefined[["alpha", "beta"].index(result["query_id"])]["query_hits"]
        }


def test_regular_and_indexed_refinement_propagate_matrix_scoring():
    database = SequenceDatabase([{"id": "gapped", "sequence": "HPEART"}])
    queries = [{"id": "query", "sequence": "HEART"}]

    regular = multi_search(
        queries,
        database,
        k=1,
        threshold=1,
        indexed=False,
        refinement=True,
        gap_penalty=-4,
        matrix="BLOSUM62",
    )
    indexed = multi_search(
        queries,
        database,
        k=1,
        threshold=1,
        indexed=True,
        refinement=True,
        gap_penalty=-4,
        matrix="BLOSUM62",
    )

    assert regular[0]["query_hits"][0]["sw_score"] == 23.0
    assert indexed[0]["query_hits"][0]["sw_score"] == 23.0


def test_indexed_matrix_refinement_ignores_simple_scores():
    database = SequenceDatabase([
        {"id": "exact", "sequence": "HEART"},
        {"id": "gapped", "sequence": "HPEART"},
    ])
    queries = [{"id": "query", "sequence": "HEART"}]

    first = multi_search(
        queries,
        database,
        k=1,
        threshold=1,
        indexed=True,
        refinement=True,
        match_score=99,
        mismatch_score=-99,
        gap_penalty=-4,
        matrix="BLOSUM62",
    )
    second = multi_search(
        queries,
        database,
        k=1,
        threshold=1,
        indexed=True,
        refinement=True,
        match_score=-50,
        mismatch_score=50,
        gap_penalty=-4,
        matrix="BLOSUM62",
    )

    assert [hit["id"] for hit in first[0]["query_hits"]] == ["exact", "gapped"]
    assert [hit["sw_score"] for hit in first[0]["query_hits"]] == [27.0, 23.0]
    assert second == first


def test_indexed_refinement_propagates_gap_penalty_for_gapped_alignment():
    database = SequenceDatabase([{"id": "gapped", "sequence": "ATTA"}])
    queries = [{"id": "query", "sequence": "ATA"}]

    stronger_penalty = multi_search(
        queries,
        database,
        k=1,
        threshold=1,
        indexed=True,
        refinement=True,
        match_score=2,
        mismatch_score=-1,
        gap_penalty=-2,
    )
    weaker_penalty = multi_search(
        queries,
        database,
        k=1,
        threshold=1,
        indexed=True,
        refinement=True,
        match_score=2,
        mismatch_score=-1,
        gap_penalty=-1,
    )

    assert stronger_penalty[0]["query_hits"][0]["sw_score"] == 4
    assert weaker_penalty[0]["query_hits"][0]["sw_score"] == 5


@pytest.mark.parametrize("indexed", [False, True])
def test_multi_search_scoring_is_inactive_without_refinement(indexed):
    database = SequenceDatabase([
        {"id": "first", "sequence": "AAAA"},
        {"id": "second", "sequence": "AATA"},
    ])
    queries = [{"id": "query", "sequence": "AAAA"}]

    default = multi_search(
        queries,
        database,
        k=1,
        threshold=1,
        indexed=indexed,
        refinement=False,
    )
    custom = multi_search(
        queries,
        database,
        k=1,
        threshold=1,
        indexed=indexed,
        refinement=False,
        match_score=100,
        mismatch_score=100,
        gap_penalty=100,
        matrix="NOT_LOADED_WHEN_REFINEMENT_IS_DISABLED",
    )

    assert custom == default
    assert all(
        "sw_score" not in hit and "best_positions" not in hit
        for hit in custom[0]["query_hits"]
    )


def test_indexed_refinement_handles_empty_queries_and_empty_hit_lists():
    database = SequenceDatabase([{"id": "hit", "sequence": "AAAA"}])

    assert multi_search(
        [],
        database,
        k=2,
        threshold=1,
        indexed=True,
        refinement=True,
    ) == []

    results = multi_search(
        [{"id": "no_hit", "sequence": "CCCC"}],
        database,
        k=2,
        threshold=1,
        indexed=True,
        refinement=True,
    )

    assert results == [{
        "query_id": "no_hit",
        "query_sequence": "CCCC",
        "query_hits": [],
    }]
