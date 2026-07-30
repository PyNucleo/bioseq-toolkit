import pytest

from bioseq.fasta_io import (
    read_fasta_records,
    read_fasta_sequences_only,
    parse_uniprot_fasta,
    parse_generic_fasta,
    parse_fasta_header,
)


PROTEIN_SEQUENCE = (
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTN"
    "AVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
)


def test_parse_uniprot_fasta_header():
    header = ">sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens OX=9606 GN=HBA1 PE=1 SV=2"

    result = parse_uniprot_fasta(header)

    assert result == {
        "id": "P69905",
        "db": "sp",
        "accession": "P69905",
        "entry_name": "HBA_HUMAN",
        "description": "Hemoglobin subunit alpha",
        "header": header,
    }


def test_parse_generic_fasta_header():
    header = ">seq1 Hemoglobin subunit alpha OS=Homo sapiens OX=9606 GN=HBA1 PE=1 SV=2"

    result = parse_generic_fasta(header)

    assert result == {
        "id": "seq1",
        "db": None,
        "accession": None,
        "entry_name": None,
        "description": "Hemoglobin subunit alpha OS=Homo sapiens OX=9606 GN=HBA1 PE=1 SV=2",
        "header": header,
    }


def test_parse_fasta_header_routes_uniprot_headers():
    header = ">sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens"

    result = parse_fasta_header(header)

    assert result["id"] == "P69905"
    assert result["db"] == "sp"
    assert result["accession"] == "P69905"
    assert result["entry_name"] == "HBA_HUMAN"
    assert result["description"] == "Hemoglobin subunit alpha"


def test_parse_fasta_header_routes_generic_headers():
    header = ">seq1 some generic sequence"

    result = parse_fasta_header(header)

    assert result["id"] == "seq1"
    assert result["db"] is None
    assert result["accession"] is None
    assert result["entry_name"] is None
    assert result["description"] == "some generic sequence"


def test_read_fasta_records_parses_uniprot_record(tmp_path):
    fasta = tmp_path / "uniprot.fasta"

    fasta.write_text(
        ">sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens OX=9606 GN=HBA1 PE=1 SV=2\n"
        f"{PROTEIN_SEQUENCE}\n"
    )

    result = read_fasta_records(fasta)[0]

    assert result == {
        "id": "P69905",
        "db": "sp",
        "accession": "P69905",
        "entry_name": "HBA_HUMAN",
        "description": "Hemoglobin subunit alpha",
        "header": ">sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens OX=9606 GN=HBA1 PE=1 SV=2",
        "sequence": PROTEIN_SEQUENCE,
    }


def test_read_fasta_records_parses_generic_record(tmp_path):
    fasta = tmp_path / "generic.fasta"

    fasta.write_text(
        ">seq1 Hemoglobin subunit alpha OS=Homo sapiens OX=9606 GN=HBA1 PE=1 SV=2\n"
        f"{PROTEIN_SEQUENCE}\n"
    )

    result = read_fasta_records(fasta)[0]

    assert result == {
        "id": "seq1",
        "db": None,
        "accession": None,
        "entry_name": None,
        "description": "Hemoglobin subunit alpha OS=Homo sapiens OX=9606 GN=HBA1 PE=1 SV=2",
        "header": ">seq1 Hemoglobin subunit alpha OS=Homo sapiens OX=9606 GN=HBA1 PE=1 SV=2",
        "sequence": PROTEIN_SEQUENCE,
    }


def test_read_fasta_records_joins_multiline_sequences(tmp_path):
    fasta = tmp_path / "multiline.fasta"

    fasta.write_text(
        ">seq1 multiline sequence\n"
        "ATGC\n"
        "GGTA\n"
        "CCAA\n"
    )

    result = read_fasta_records(fasta)[0]

    assert result["id"] == "seq1"
    assert result["sequence"] == "ATGCGGTACCAA"


def test_read_fasta_sequences_only_returns_sequences(tmp_path):
    fasta = tmp_path / "sequences_only.fasta"

    fasta.write_text(
        ">seq1 first sequence\n"
        "ATGC\n"
        ">seq2 second sequence\n"
        "GGCC\n"
    )

    result = read_fasta_sequences_only(fasta)

    assert result == ["ATGC", "GGCC"]


def test_read_fasta_records_ignores_blank_lines_and_joins_sequences(tmp_path):
    fasta = tmp_path / "blank_lines.fasta"

    fasta.write_text(
        "\n"
        "   \n"
        ">seq1 first sequence\n"
        "\n"
        "ATGC\n"
        "\t\n"
        "GGTA\n"
        "\n"
        ">seq2 second sequence\n"
        "CCAA\n"
        "TTGG\n"
        "   \n"
    )

    result = read_fasta_records(fasta)

    assert [(record["id"], record["description"], record["sequence"]) for record in result] == [
        ("seq1", "first sequence", "ATGCGGTA"),
        ("seq2", "second sequence", "CCAATTGG"),
    ]


@pytest.mark.parametrize("header", [">\n", ">   \n", ">>\n"])
def test_read_fasta_records_rejects_empty_header_ids(tmp_path, header):
    fasta = tmp_path / "empty_header.fasta"
    fasta.write_text(f"\n{header}ATGC\n")

    with pytest.raises(ValueError) as error:
        read_fasta_records(fasta)

    message = str(error.value).lower()
    assert "line 2" in message
    assert "record id" in message
    assert "empty" in message


def test_read_fasta_records_rejects_sequence_before_first_header(tmp_path):
    fasta = tmp_path / "sequence_before_header.fasta"
    fasta.write_text("\n\t\nATGC\n>seq1\nGGTA\n")

    with pytest.raises(ValueError) as error:
        read_fasta_records(fasta)

    message = str(error.value).lower()
    assert "line 3" in message
    assert "sequence" in message
    assert "before the first header" in message


@pytest.mark.parametrize("sequence", ["AT GC", "AT\tGC", " ATGC", "ATGC "])
def test_read_fasta_records_rejects_whitespace_in_sequence_lines(tmp_path, sequence):
    fasta = tmp_path / "sequence_whitespace.fasta"
    fasta.write_text(f">seq1\n{sequence}\n")

    with pytest.raises(ValueError) as error:
        read_fasta_records(fasta)

    message = str(error.value).lower()
    assert "line 2" in message
    assert "sequence" in message
    assert "whitespace" in message


def test_read_fasta_records_rejects_empty_sequence_before_next_header(tmp_path):
    fasta = tmp_path / "empty_sequence_before_next_header.fasta"
    fasta.write_text(
        ">seq1\n"
        "\n"
        "\n"
        ">seq2\n"
        "ATGC\n"
    )

    with pytest.raises(ValueError) as error:
        read_fasta_records(fasta)

    message = str(error.value).lower()
    assert "line 1" in message
    assert "sequence" in message
    assert "empty" in message


def test_read_fasta_records_rejects_empty_final_sequence(tmp_path):
    fasta = tmp_path / "empty_final_sequence.fasta"
    fasta.write_text(
        ">seq1\n"
        "\n"
        "\n"
    )

    with pytest.raises(ValueError) as error:
        read_fasta_records(fasta)

    message = str(error.value).lower()
    assert "line 1" in message
    assert "sequence" in message
    assert "empty" in message


@pytest.mark.parametrize("contents", ["", "\n \n\t\n"])
def test_read_fasta_records_returns_no_records_for_empty_or_blank_files(tmp_path, contents):
    fasta = tmp_path / "blank.fasta"
    fasta.write_text(contents)

    assert read_fasta_records(fasta) == []


def test_read_fasta_records_accepts_crlf_input(tmp_path):
    fasta = tmp_path / "crlf.fasta"
    fasta.write_bytes(b">seq1 description\r\nATGC\r\nGGTA\r\n")

    result = read_fasta_records(fasta)

    assert result[0]["header"] == ">seq1 description"
    assert result[0]["sequence"] == "ATGCGGTA"
