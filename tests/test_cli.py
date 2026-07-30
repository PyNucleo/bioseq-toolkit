import json
import subprocess
import pytest
import bioseq.cli as cli
import sys


def run_cli(*args):
    """
    Run the bioseq CLI as a real subprocess.

    This tests the command-line interface the same way a user would run it
    from the terminal.
    """
    return subprocess.run(
        [sys.executable, "-m", "bioseq.cli", *args],
        capture_output=True,
        text=True,
    )


def test_cli_main_help_runs():
    result = run_cli("--help")

    assert result.returncode == 0
    assert "search" in result.stdout
    assert "align-local" in result.stdout
    assert "align-global" in result.stdout


def test_cli_search_help_runs():
    result = run_cli("search", "--help")

    assert result.returncode == 0
    assert "--query" in result.stdout
    assert "--database" in result.stdout
    assert "--kmer-size" in result.stdout
    assert "--threshold" in result.stdout
    assert "--refine" in result.stdout


def test_cli_align_local_help_runs():
    result = run_cli("align-local", "--help")

    assert result.returncode == 0
    assert "--sequence1" in result.stdout
    assert "--sequence2" in result.stdout
    assert "--match" in result.stdout
    assert "--mismatch" in result.stdout
    assert "--gap-penalty" in result.stdout
    assert "--matrix" in result.stdout


def test_cli_align_global_help_runs():
    result = run_cli("align-global", "--help")

    assert result.returncode == 0
    assert "--sequence1" in result.stdout
    assert "--sequence2" in result.stdout
    assert "--match" in result.stdout
    assert "--mismatch" in result.stdout
    assert "--gap-penalty" in result.stdout
    assert "--matrix" in result.stdout


def test_cli_search_returns_kmer_hits(tmp_path):
    fasta = tmp_path / "test.fasta"

    fasta.write_text(
        ">seq1 best hit\n"
        "ATGCGT\n"
        ">seq2 also good\n"
        "ATGCGA\n"
        ">seq3 weak hit\n"
        "ATGAAA\n"
        ">seq4 no hit\n"
        "GGGGGG\n"
    )

    result = run_cli(
        "search",
        "-q", "ATGCG",
        "-d", str(fasta),
        "-k", "3",
        "-t", "2",
        "-n", "10",
    )

    assert result.returncode == 0, result.stderr

    data = json.loads(result.stdout)

    assert isinstance(data, list)
    assert len(data) == 2

    assert data[0]["id"] == "seq1"
    assert data[0]["sequence"] == "ATGCGT"
    assert data[0]["shared_kmers"] == 3

    assert data[1]["id"] == "seq2"
    assert data[1]["sequence"] == "ATGCGA"
    assert data[1]["shared_kmers"] == 3

    assert all("sw_score" not in hit for hit in data)


def test_cli_search_with_refinement_adds_sw_scores(tmp_path):
    fasta = tmp_path / "test.fasta"

    fasta.write_text(
        ">seq1 best hit\n"
        "ATGCGT\n"
        ">seq2 also good\n"
        "ATGCGA\n"
        ">seq3 weak hit\n"
        "ATGAAA\n"
        ">seq4 no hit\n"
        "GGGGGG\n"
    )

    result = run_cli(
        "search",
        "-q", "ATGCG",
        "-d", str(fasta),
        "-k", "3",
        "-t", "1",
        "-n", "3",
        "-r",
    )

    assert result.returncode == 0, result.stderr

    data = json.loads(result.stdout)

    assert isinstance(data, list)
    assert len(data) == 3

    assert all("id" in hit for hit in data)
    assert all("sequence" in hit for hit in data)
    assert all("shared_kmers" in hit for hit in data)
    assert all("sw_score" in hit for hit in data)
    assert all("best_positions" in hit for hit in data)

    assert data == sorted(
        data,
        key=lambda hit: hit["sw_score"],
        reverse=True,
    )


def test_cli_align_local_returns_smith_waterman_result():
    result = run_cli(
        "align-local",
        "-s1", "HEART",
        "-s2", "HPEART",
        "--matrix", "BLOSUM62",
        "-g", "-4",
    )

    assert result.returncode == 0, result.stderr

    data = json.loads(result.stdout)

    assert data["algorithm"] == "Smith-Waterman"
    assert data["mode"] == "local"
    assert data["sequence_1"] == "HEART"
    assert data["sequence_2"] == "HPEART"
    assert data["score"] == 23.0
    assert data["scoring"] == {
        "gap_penalty": -4,
        "matrix": "BLOSUM62",
        "gap_model": "linear",
    }

    assert data["best_positions"] == [[5, 6]]
    assert data["num_alignments"] == 1

    alignment = data["alignments"][0]

    assert alignment["aligned_sequence_1"] == "H-EART"
    assert alignment["aligned_sequence_2"] == "HPEART"
    assert alignment["matches"] == 5
    assert alignment["gaps"] == 1
    assert alignment["identity"] == 0.833


def test_cli_align_global_returns_needleman_wunsch_result():
    result = run_cli(
        "align-global",
        "-s1", "ATGCG",
        "-s2", "ATCGA",
        "-m", "1",
        "--mismatch", "-1",
        "-g", "-2",
    )

    assert result.returncode == 0, result.stderr

    data = json.loads(result.stdout)

    assert data["algorithm"] == "Needleman-Wunsch"
    assert data["mode"] == "global"
    assert data["sequence_1"] == "ATGCG"
    assert data["sequence_2"] == "ATCGA"
    assert data["score"] == 0
    assert data["scoring"] == {
        "match": 1,
        "mismatch": -1,
        "gap_penalty": -2,
        "matrix": None,
        "gap_model": "linear",
    }

    assert data["num_alignments"] == 1

    alignment = data["alignments"][0]

    assert alignment["aligned_sequence_1"] == "ATGCG-"
    assert alignment["aligned_sequence_2"] == "AT-CGA"
    assert alignment["matches"] == 4
    assert alignment["gaps"] == 2
    assert alignment["identity"] == 0.667


def test_cli_rejects_missing_required_search_arguments():
    result = run_cli(
        "search",
        "-q", "ATGCG",
    )

    assert result.returncode != 0
    assert "database" in result.stderr.lower()


def test_cli_rejects_unknown_command():
    result = run_cli("not-a-command")

    assert result.returncode != 0
    assert "invalid choice" in result.stderr.lower()

@pytest.mark.parametrize(
    ("extra_args", "expected_full_header"),
    [
        ([], False),
        (["--full-header"], True),
    ],
)
def test_cli_fetch_uniprot_passes_full_header_flag(
    monkeypatch,
    tmp_path,
    extra_args,
    expected_full_header,
):
    input_path = tmp_path / "accessions.txt"
    output_path = tmp_path / "output.fasta"

    fetched_result = {
        "records": [
            {
                "id": "P12345",
                "accession": "P12345",
                "header": ">sp|P12345|TEST_ENTRY Example protein",
                "sequence": "ATGC",
            }
        ],
        "failed": [],
    }

    captured = {}

    def fake_fetch_uniprot_sequences(file_path, strict):
        captured["fetch_args"] = (file_path, strict)
        return fetched_result

    def fake_write_fasta_records(records, result_path, full_header):
        captured["write_args"] = (records, result_path, full_header)

    monkeypatch.setattr(
        cli,
        "fetch_uniprot_sequences",
        fake_fetch_uniprot_sequences,
    )
    monkeypatch.setattr(
        cli,
        "write_fasta_records",
        fake_write_fasta_records,
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "bioseq",
            "fetch-uniprot",
            "-f",
            str(input_path),
            "-o",
            str(output_path),
            *extra_args,
        ],
    )

    cli.main()

    assert captured["fetch_args"] == (str(input_path), False)
    assert captured["write_args"] == (
        fetched_result["records"],
        str(output_path),
        expected_full_header,
    )