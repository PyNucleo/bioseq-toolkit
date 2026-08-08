import json

import pytest

from bioseq.pipelines.translation_pipeline import process_fasta_sequences


def test_process_fasta_sequences_accepts_mixed_case_dna_and_preserves_outputs(tmp_path):
    fasta = tmp_path / "valid.fasta"
    fasta.write_text(
        ">seq1 first record\n"
        "atgc\n"
        ">seq2\n"
        "ATGGCCATT\n"
    )

    result = process_fasta_sequences(fasta)

    assert result == {
        "accepted": [
            {
                "id": "seq1",
                "record_position": 1,
                "sequence": "ATGC",
                "length": 4,
                "gc_content": 50.0,
                "transcribed_strand": "UACG",
                "amino_acid_chain": "Y",
            },
            {
                "id": "seq2",
                "record_position": 2,
                "sequence": "ATGGCCATT",
                "length": 9,
                "gc_content": 44.44,
                "transcribed_strand": "UACCGGUAA",
                "amino_acid_chain": "YR",
            },
        ],
        "rejected": [],
        "summary": {
            "total_records": 2,
            "accepted_records": 2,
            "rejected_records": 0,
        },
    }
    assert isinstance(result["accepted"][0]["amino_acid_chain"], str)
    json.dumps(result)


def test_process_fasta_sequences_reports_interleaved_invalid_records(tmp_path):
    fasta = tmp_path / "mixed.fasta"
    fasta.write_text(
        ">bad-first\n"
        "AT-N\n"
        ">good-second\n"
        "ATGC\n"
        ">bad-third\n"
        "GG-R\n"
        ">good-fourth\n"
        "CCAA\n"
    )

    result = process_fasta_sequences(fasta)

    assert [record["record_position"] for record in result["accepted"]] == [2, 4]
    assert [record["record_position"] for record in result["rejected"]] == [1, 3]
    assert result["rejected"][0] == {
        "id": "bad-first",
        "record_position": 1,
        "sequence": "AT-N",
        "reason_code": "unsupported_dna_symbols",
        "reason": "Sequence contains unsupported DNA symbols.",
        "invalid_positions": [3, 4],
        "invalid_symbols": ["-", "N"],
    }
    _assert_summary_invariants(result)


def test_process_fasta_sequences_reports_repeated_invalid_symbol_once(tmp_path):
    fasta = tmp_path / "repeated.fasta"
    fasta.write_text(">seq1\nA--TG-\n")

    rejected = process_fasta_sequences(fasta)["rejected"]

    assert rejected[0]["invalid_positions"] == [2, 3, 6]
    assert rejected[0]["invalid_symbols"] == ["-"]
    assert isinstance(rejected[0]["invalid_symbols"], list)


def test_process_fasta_sequences_reports_multiple_symbols_deterministically(tmp_path):
    fasta = tmp_path / "multiple.fasta"
    fasta.write_text(">seq1\nAT-NRY\n")

    rejected = process_fasta_sequences(fasta)["rejected"]

    assert rejected[0]["invalid_positions"] == [3, 4, 5, 6]
    assert rejected[0]["invalid_symbols"] == ["-", "N", "R", "Y"]


def test_process_fasta_sequences_rejects_n_and_all_invalid_records(tmp_path):
    fasta = tmp_path / "all-invalid.fasta"
    fasta.write_text(">n-record\nATN\n>gap-record\n--\n")

    result = process_fasta_sequences(fasta)

    assert result["accepted"] == []
    assert [record["id"] for record in result["rejected"]] == ["n-record", "gap-record"]
    _assert_summary_invariants(result)


def test_process_fasta_sequences_returns_empty_result_for_empty_fasta(tmp_path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("")

    assert process_fasta_sequences(fasta) == {
        "accepted": [],
        "rejected": [],
        "summary": {
            "total_records": 0,
            "accepted_records": 0,
            "rejected_records": 0,
        },
    }


def test_process_fasta_sequences_propagates_malformed_empty_record_error(tmp_path):
    fasta = tmp_path / "malformed.fasta"
    fasta.write_text(">seq1\n>seq2\nATGC\n")

    with pytest.raises(ValueError, match="sequence is empty"):
        process_fasta_sequences(fasta)


def test_process_fasta_sequences_keeps_multiline_sequence_as_one_record(tmp_path):
    fasta = tmp_path / "multiline.fasta"
    fasta.write_text(">seq1\natgc\nggta\n")

    accepted = process_fasta_sequences(fasta)["accepted"]

    assert accepted[0]["record_position"] == 1
    assert accepted[0]["sequence"] == "ATGCGGTA"


def _assert_summary_invariants(result):
    summary = result["summary"]
    assert summary["accepted_records"] == len(result["accepted"])
    assert summary["rejected_records"] == len(result["rejected"])
    assert summary["total_records"] == (
        summary["accepted_records"] + summary["rejected_records"]
    )
