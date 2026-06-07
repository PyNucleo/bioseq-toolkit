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