import pytest

from bioseq.fasta_io import write_fasta_records


def test_write_fasta_records_with_full_headers(tmp_path):
    records = [
        {
            "id": "Q9HUU8",
            "accession": "Q9HUU8",
            "header": ">sp|Q9HUU8|URE3_PSEAE Urease subunit gamma OS=Pseudomonas aeruginosa",
            "sequence": "MDLSPREKDKLLIFTAGLLAERRLARGLKLNYPEA",
        },
        {
            "id": "P02748",
            "accession": "P02748",
            "header": ">sp|P02748|CO9_HUMAN Complement component C9 OS=Homo sapiens",
            "sequence": "MSACRSFAVAICILEISILTAQYTTSYDPELTESSG",
        },
    ]

    output_file = tmp_path / "records.fasta"

    write_fasta_records(records, output_file, full_header=True)

    result = output_file.read_text()

    expected = (
        ">sp|Q9HUU8|URE3_PSEAE Urease subunit gamma OS=Pseudomonas aeruginosa\n"
        "MDLSPREKDKLLIFTAGLLAERRLARGLKLNYPEA\n"
        ">sp|P02748|CO9_HUMAN Complement component C9 OS=Homo sapiens\n"
        "MSACRSFAVAICILEISILTAQYTTSYDPELTESSG\n"
    )

    assert result == expected


def test_write_fasta_records_with_simple_headers_uses_accession(tmp_path):
    records = [
        {
            "id": "Q9HUU8",
            "accession": "Q9HUU8",
            "header": ">sp|Q9HUU8|URE3_PSEAE Urease subunit gamma OS=Pseudomonas aeruginosa",
            "sequence": "MDLSPREKDKLLIFTAGLLAERRLARGLKLNYPEA",
        },
        {
            "id": "P02748",
            "accession": "P02748",
            "header": ">sp|P02748|CO9_HUMAN Complement component C9 OS=Homo sapiens",
            "sequence": "MSACRSFAVAICILEISILTAQYTTSYDPELTESSG",
        },
    ]

    output_file = tmp_path / "records_simple.fasta"

    write_fasta_records(records, output_file, full_header=False)

    result = output_file.read_text()

    expected = (
        ">Q9HUU8\n"
        "MDLSPREKDKLLIFTAGLLAERRLARGLKLNYPEA\n"
        ">P02748\n"
        "MSACRSFAVAICILEISILTAQYTTSYDPELTESSG\n"
    )

    assert result == expected


def test_write_fasta_records_simple_header_falls_back_to_id(tmp_path):
    records = [
        {
            "id": "seq1",
            "header": ">seq1 example sequence",
            "sequence": "ATGCGT",
        }
    ]

    output_file = tmp_path / "generic.fasta"

    write_fasta_records(records, output_file, full_header=False)

    result = output_file.read_text()

    expected = (
        ">seq1\n"
        "ATGCGT\n"
    )

    assert result == expected


def test_write_fasta_records_does_not_duplicate_header_marker(tmp_path):
    records = [
        {
            "id": "Q9HUU8",
            "accession": "Q9HUU8",
            "header": ">sp|Q9HUU8|URE3_PSEAE Urease subunit gamma",
            "sequence": "MDLSPREK",
        }
    ]

    output_file = tmp_path / "records.fasta"

    write_fasta_records(records, output_file, full_header=True)

    result = output_file.read_text()

    assert result.startswith(">sp|Q9HUU8")
    assert not result.startswith(">>")


def test_write_fasta_records_short_header_prefers_different_accession(tmp_path):
    output_file = tmp_path / "accession.fasta"

    write_fasta_records(
        [{"id": "ENTRY_NAME", "accession": "P12345", "sequence": "ATGC"}],
        output_file,
        full_header=False,
    )

    assert output_file.read_text() == ">P12345\nATGC\n"


@pytest.mark.parametrize(
    "record",
    [
        {"id": "seq1", "sequence": "ATGC"},
        {"id": "seq1", "accession": None, "sequence": "ATGC"},
        {"id": "seq1", "accession": "", "sequence": "ATGC"},
    ],
)
def test_write_fasta_records_short_header_falls_back_to_id(tmp_path, record):
    output_file = tmp_path / "generic.fasta"

    write_fasta_records([record], output_file, full_header=False)

    assert output_file.read_text() == ">seq1\nATGC\n"


def test_write_fasta_records_selects_identifiers_per_record(tmp_path):
    output_file = tmp_path / "mixed.fasta"
    records = [
        {"id": "ENTRY_A", "accession": "P12345", "sequence": "ATGC"},
        {"id": "seq2", "accession": None, "sequence": "GGCC"},
        {"id": "seq3", "sequence": "TTAA"},
    ]

    write_fasta_records(records, output_file, full_header=False)

    assert output_file.read_text() == ">P12345\nATGC\n>seq2\nGGCC\n>seq3\nTTAA\n"


def test_write_fasta_records_short_header_removes_existing_marker(tmp_path):
    output_file = tmp_path / "marker.fasta"

    write_fasta_records(
        [{"id": "seq1", "accession": ">P12345", "sequence": "ATGC"}],
        output_file,
        full_header=False,
    )

    assert output_file.read_text() == ">P12345\nATGC\n"


def test_write_fasta_records_rejects_missing_short_identifier(tmp_path):
    with pytest.raises(ValueError, match=r"Record 1.*usable accession or ID"):
        write_fasta_records(
            [{"id": "", "accession": None, "sequence": "ATGC"}],
            tmp_path / "missing-identifier.fasta",
            full_header=False,
        )


def test_write_fasta_records_rejects_missing_sequence_with_record_number(tmp_path):
    records = [
        {"id": "seq1", "sequence": "ATGC"},
        {"id": "seq2"},
    ]

    with pytest.raises(ValueError, match=r"Record 2.*sequence"):
        write_fasta_records(records, tmp_path / "missing-sequence.fasta", full_header=False)


def test_write_fasta_records_rejects_missing_full_header(tmp_path):
    with pytest.raises(ValueError, match=r"Record 1.*header"):
        write_fasta_records(
            [{"id": "seq1", "sequence": "ATGC"}],
            tmp_path / "missing-header.fasta",
            full_header=True,
        )
