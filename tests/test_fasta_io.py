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