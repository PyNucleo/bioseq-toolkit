from types import SimpleNamespace
from unittest.mock import Mock, call

import pytest
import requests

import bioseq.fasta_io as fasta_io


def make_response(status_code, text=""):
    return SimpleNamespace(status_code=status_code, text=text)


def make_uniprot_fasta(accession, entry_name):
    return (
        f">sp|{accession}|{entry_name} Example protein OS=Homo sapiens\n"
        "MABCDE\n"
    )


def test_fetch_uniprot_success_matches_local_parser_and_returns_flat_record(
    monkeypatch,
    tmp_path,
):
    fasta_body = make_uniprot_fasta("P12345", "TEST_HUMAN")
    local_fasta = tmp_path / "record.fasta"
    local_fasta.write_text(fasta_body)
    accession_file = tmp_path / "accessions.txt"
    accession_file.write_text("P12345\n")
    mock_get = Mock(return_value=make_response(200, fasta_body))
    monkeypatch.setattr(fasta_io.requests, "get", mock_get)

    local_record = fasta_io.read_fasta_records(local_fasta)[0]
    result = fasta_io.fetch_uniprot_sequences(accession_file, strict=False)

    assert result["failed"] == []
    assert len(result["records"]) == 1
    assert isinstance(result["records"][0], dict)
    assert result["records"][0] == local_record
    assert result["records"][0] == {
        "id": "P12345",
        "db": "sp",
        "accession": "P12345",
        "entry_name": "TEST_HUMAN",
        "description": "Example protein",
        "header": ">sp|P12345|TEST_HUMAN Example protein OS=Homo sapiens",
        "sequence": "MABCDE",
    }
    mock_get.assert_called_once_with(
        "https://www.uniprot.org/uniprotkb/P12345.fasta"
    )


def test_fetch_uniprot_multiple_accessions_skip_blank_lines_and_preserve_order(
    monkeypatch,
    tmp_path,
):
    accession_file = tmp_path / "accessions.txt"
    accession_file.write_text("P11111\n\nP22222\n")
    mock_get = Mock(
        side_effect=[
            make_response(200, make_uniprot_fasta("P11111", "FIRST_HUMAN")),
            make_response(200, make_uniprot_fasta("P22222", "SECOND_HUMAN")),
        ]
    )
    monkeypatch.setattr(fasta_io.requests, "get", mock_get)

    result = fasta_io.fetch_uniprot_sequences(accession_file, strict=False)

    assert result["failed"] == []
    assert len(result["records"]) == 2
    assert all(isinstance(record, dict) for record in result["records"])
    assert [record["accession"] for record in result["records"]] == [
        "P11111",
        "P22222",
    ]
    assert mock_get.call_args_list == [
        call("https://www.uniprot.org/uniprotkb/P11111.fasta"),
        call("https://www.uniprot.org/uniprotkb/P22222.fasta"),
    ]


def test_fetch_uniprot_non_strict_continues_after_http_failure(
    monkeypatch,
    tmp_path,
):
    accession_file = tmp_path / "accessions.txt"
    accession_file.write_text("BAD\nGOOD\n")
    mock_get = Mock(
        side_effect=[
            make_response(404),
            make_response(200, make_uniprot_fasta("GOOD", "GOOD_HUMAN")),
        ]
    )
    monkeypatch.setattr(fasta_io.requests, "get", mock_get)

    result = fasta_io.fetch_uniprot_sequences(accession_file, strict=False)

    assert len(result["failed"]) == 1
    failure = result["failed"][0]
    assert failure["accession"] == "BAD"
    assert failure["status_code"] == 404
    assert failure["reason"]
    assert [record["accession"] for record in result["records"]] == ["GOOD"]
    assert mock_get.call_count == 2


def test_fetch_uniprot_strict_stops_after_http_failure(monkeypatch, tmp_path):
    accession_file = tmp_path / "accessions.txt"
    accession_file.write_text("BAD\nGOOD\n")
    mock_get = Mock(
        side_effect=[
            make_response(404),
            make_response(200, make_uniprot_fasta("GOOD", "GOOD_HUMAN")),
        ]
    )
    monkeypatch.setattr(fasta_io.requests, "get", mock_get)

    with pytest.raises(ValueError) as error:
        fasta_io.fetch_uniprot_sequences(accession_file, strict=True)

    assert "BAD" in str(error.value)
    assert "404" in str(error.value)
    assert mock_get.call_count == 1


def test_fetch_uniprot_non_strict_continues_after_network_error(
    monkeypatch,
    tmp_path,
):
    accession_file = tmp_path / "accessions.txt"
    accession_file.write_text("BAD\nGOOD\n")
    mock_get = Mock(
        side_effect=[
            requests.ConnectionError("connection lost"),
            make_response(200, make_uniprot_fasta("GOOD", "GOOD_HUMAN")),
        ]
    )
    monkeypatch.setattr(fasta_io.requests, "get", mock_get)

    result = fasta_io.fetch_uniprot_sequences(accession_file, strict=False)

    assert len(result["failed"]) == 1
    failure = result["failed"][0]
    assert failure["accession"] == "BAD"
    assert failure["status_code"] is None
    assert "connection lost" in failure["reason"]
    assert [record["accession"] for record in result["records"]] == ["GOOD"]
    assert mock_get.call_count == 2


def test_fetch_uniprot_strict_stops_after_network_error(monkeypatch, tmp_path):
    accession_file = tmp_path / "accessions.txt"
    accession_file.write_text("BAD\nGOOD\n")
    network_error = requests.Timeout("timed out")
    mock_get = Mock(
        side_effect=[
            network_error,
            make_response(200, make_uniprot_fasta("GOOD", "GOOD_HUMAN")),
        ]
    )
    monkeypatch.setattr(fasta_io.requests, "get", mock_get)

    with pytest.raises(ValueError) as error:
        fasta_io.fetch_uniprot_sequences(accession_file, strict=True)

    assert "BAD" in str(error.value)
    assert "timed out" in str(error.value)
    assert error.value.__cause__ is network_error
    assert mock_get.call_count == 1


def test_fetch_uniprot_non_strict_continues_after_empty_http_200(
    monkeypatch,
    tmp_path,
):
    accession_file = tmp_path / "accessions.txt"
    accession_file.write_text("EMPTY\nGOOD\n")
    mock_get = Mock(
        side_effect=[
            make_response(200, ""),
            make_response(200, make_uniprot_fasta("GOOD", "GOOD_HUMAN")),
        ]
    )
    monkeypatch.setattr(fasta_io.requests, "get", mock_get)

    result = fasta_io.fetch_uniprot_sequences(accession_file, strict=False)

    assert len(result["failed"]) == 1
    failure = result["failed"][0]
    assert failure["accession"] == "EMPTY"
    assert failure["status_code"] == 200
    assert failure["reason"]
    assert [record["accession"] for record in result["records"]] == ["GOOD"]
    assert mock_get.call_count == 2


def test_fetch_uniprot_strict_stops_after_empty_http_200(monkeypatch, tmp_path):
    accession_file = tmp_path / "accessions.txt"
    accession_file.write_text("EMPTY\nGOOD\n")
    mock_get = Mock(
        side_effect=[
            make_response(200, ""),
            make_response(200, make_uniprot_fasta("GOOD", "GOOD_HUMAN")),
        ]
    )
    monkeypatch.setattr(fasta_io.requests, "get", mock_get)

    with pytest.raises(ValueError, match="Empty or invalid FASTA response"):
        fasta_io.fetch_uniprot_sequences(accession_file, strict=True)

    assert mock_get.call_count == 1
