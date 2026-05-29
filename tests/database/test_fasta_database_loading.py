# tests/database/test_fasta_database_loading.py

from bioseq.pipelines.search_pipeline import search


def test_search_accepts_fasta_file(tmp_path):

    fasta = tmp_path / "test.fasta"

    fasta.write_text(
        ">seqA\n"
        "ATGCGT\n"
        ">seqB\n"
        "GGGGGG\n"
    )

    results = search(
        query="ATGCG",
        database=str(fasta),
        k=3,
        threshold=1
    )

    assert len(results) > 0
    assert results[0]["sequence"] == "ATGCGT"
    assert "id" in results[0]