import requests

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

def parse_header_sequence_from_string(text):

    text = text.strip()

    split_text = text.splitlines()

    header_seq_dict = {
        "header": split_text[0],
        "sequence": "".join(split_text[1:])
    }

    return header_seq_dict

def fetch_uniprot_sequences(accession_file, strict=False):

    parsed_sequences = []
    failed_accessions = []

    with open(accession_file, "r") as f:

        for accession in f:

            accession = accession.strip()

            if not accession:
                continue

            temp_url = f"https://www.uniprot.org/uniprotkb/{accession}.fasta"

            response = requests.get(temp_url)

            if response.status_code == 200:
                fasta_data = response.text
                header_and_sequence = parse_header_sequence_from_string(fasta_data)

                header = header_and_sequence["header"]
                sequence = header_and_sequence["sequence"]

                entry_data = {
                    **parse_fasta_header(header),
                    "sequence": sequence
                }

                parsed_sequences.append(entry_data)

            else:
                if strict:
                    raise ValueError(
                        f"Failed to fetch sequence for accession {accession}. "
                        f"Status code: {response.status_code}"
                    )
                
                failed_accessions.append({
                    "accession": accession,
                    "status_code": response.status_code
                })

    return {
        "records": parsed_sequences, 
        "failed": failed_accessions
    }

def write_fasta_records(records, output_path, full_header):

    if full_header:
        word = "header"
        first_char = ""
    else:
        first_char = ">"

        if "accession" in records:
            word = "accession"
        else:
            word = "id"
        
    with open(output_path, "w", newline="") as file:

        for record in records:
            file.write(first_char + record[word] + "\n")
            file.write(record["sequence"] + "\n")