def read_fasta_records(file):
    FASTA = open(file, "r")
    records = []
    record = None
    current_seq = ""

    for line in FASTA:

        line = line.rstrip()

        if line.startswith(">"):

            if record is not None:
                record["sequence"] = current_seq
                records.append(record)

            record = parse_fasta_header(line)
            record["sequence"] = ""
            

            current_seq = ""

        else:

            current_seq += line

    if record is not None:

        record["sequence"] = current_seq
        records.append(record)

    return records

def parse_fasta_header(header):
    clean_header = header.lstrip(">")

    if clean_header.startswith(("sp|", "tr|")):
        return parse_uniprot_fasta(header)

    return parse_generic_fasta(header)



def read_fasta_sequences_only(FILE):

    records = read_fasta_records(FILE)

    sequences = [record["sequence"] for record in records]

    return sequences

def parse_uniprot_fasta(header):
    """
    Parse a UniProt-style FASTA header.

    Example:
    >sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens OX=9606

    Returns:
    {
        "id": "P69905",
        "db": "sp",
        "accession": "P69905",
        "entry_name": "HBA_HUMAN",
        "description": "Hemoglobin subunit alpha",
        "header": original_header
    }
    """
    clean_header = header.lstrip(">")

    parts = clean_header.split("|", 2)

    if len(parts) != 3:
        return parse_generic_fasta(header)

    db = parts[0]
    accession = parts[1]
    remainder = parts[2]

    tokens = remainder.split()

    if not tokens:
        return {
            "id": accession,
            "db": db,
            "accession": accession,
            "entry_name": "",
            "description": "",
            "header": header
        }

    entry_name = tokens[0]
    description_tokens = []

    for token in tokens[1:]:
        if token.startswith("OS="):
            break

        description_tokens.append(token)

    description = " ".join(description_tokens)

    return {
        "id": accession,
        "db": db,
        "accession": accession,
        "entry_name": entry_name,
        "description": description,
        "header": header
    }


def parse_generic_fasta(header):
    """
    Parse a generic FASTA header.

    Example:
    >seqA some description

    Returns:
    {
        "id": "seqA",
        "db": None,
        "accession": None,
        "entry_name": None,
        "description": "some description",
        "header": original_header
    }
    """
    clean_header = header.lstrip(">")
    tokens = clean_header.split()

    if not tokens:
        return {
            "id": "",
            "db": None,
            "accession": None,
            "entry_name": None,
            "description": "",
            "header": header
        }

    record_id = tokens[0]
    description = " ".join(tokens[1:])

    return {
        "id": record_id,
        "db": None,
        "accession": None,
        "entry_name": None,
        "description": description,
        "header": header
    }


