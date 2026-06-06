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