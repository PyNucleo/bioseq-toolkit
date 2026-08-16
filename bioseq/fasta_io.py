import requests

def _parse_fasta_records(fasta_data):

    records = []
    record = None
    current_seq = ""
    record_starting_line = None

    for line_number, line in enumerate(fasta_data, start=1):

        line = line.rstrip("\r\n")

        if not line.strip():
            continue

        if line.startswith(">"):
            if record is not None:
                if not current_seq:
                    raise ValueError(f'Invalid FASTA detected on line {record_starting_line}: '
                                        'sequence is empty.')
                record["sequence"] = current_seq
                records.append(record)

            record_starting_line = line_number
            record = parse_fasta_header(line)
            validate_fasta_header(record, line_number)
            record["sequence"] = ""
            

            current_seq = ""

        else:
            if record is None: #Did not find a previous header
                raise ValueError(f'Invalid FASTA detected on line {line_number}: '
                                    'Sequence data appeared before the first header.')
            if any(character.isspace() for character in line):
                raise ValueError(
                        f"Invalid FASTA sequence at line {line_number}: "
                        "sequence data contains whitespace."
                )
            current_seq += line

    if record is not None:
        if not current_seq:
            raise ValueError(f'Invalid FASTA detected on line {record_starting_line}: '
                                'sequence is empty.')
        record["sequence"] = current_seq
        records.append(record)

    return records
    
def validate_fasta_header(record, line_number=None):
    if not record.get("id"):
        location = (
            f" at line {line_number}"
            if line_number is not None
            else ""
        )

        raise ValueError(
            f"Invalid FASTA header{location}: the record ID is empty."
        )

def read_fasta_records(file_path):

    with open(file_path) as fasta_file:
        return _parse_fasta_records(fasta_file)

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

            try:
                response = requests.get(temp_url)
            except requests.RequestException as error:
                if strict:
                    raise ValueError(
                        f"Failed to fetch sequence for accession {accession}: {error}"
                    ) from error

                failed_accessions.append({
                    "accession": accession,
                    "status_code": None,
                    "reason": f"Network request failed: {error}"
                })
                continue


            if response.status_code == 200:
                fasta_data = response.text.strip()

                if not fasta_data:

                    if strict:
                        raise ValueError("Empty or invalid FASTA response. The accession may be secondary, "
                                                    "merged, demerged, obsolete, ambiguous, or otherwise not directly "
                                                    "resolvable to a current UniProt FASTA record. Tip: search this accession "
                                                    "manually on UniProt to check its entry history and possible current "
                                                    "primary accession(s).")
                    failed_accessions.append(
                        {"accession": accession,
                         "status_code": 200,
                         "reason": (
                            "Empty or invalid FASTA response. The accession may be secondary, "
                            "merged, demerged, obsolete, ambiguous, or otherwise not directly "
                            "resolvable to a current UniProt FASTA record. Tip: search this accession "
                            "manually on UniProt to check its entry history and possible current "
                            "primary accession(s)."
                         )
                        })
                    continue

                entry_data = _parse_fasta_records(fasta_data.splitlines())
                parsed_sequences.extend(entry_data)

            else:
                if strict:
                    raise ValueError(
                        f"Failed to fetch sequence for accession {accession}. "
                        f"Status code: {response.status_code}"
                    )

                failed_accessions.append({
                    "accession": accession,
                    "status_code": response.status_code,
                    "reason": "Reason: UniProt did not return a FASTA record for this accession. The accession may be invalid, retired, merged, obsolete, or incorrectly typed."
                })

    return {
        "records": parsed_sequences, 
        "failed": failed_accessions
    }

def write_fasta_records(records, output_path, full_header):

    with open(output_path, "w", newline="") as file:
        for record_number, record in enumerate(records, start=1):
            if "sequence" not in record:
                raise ValueError(f"Record {record_number} is missing required sequence field.")

            if full_header:
                header = record.get("header")
                if not isinstance(header, str) or not header:
                    raise ValueError(
                        f"Record {record_number} has a missing or unusable header."
                    )
            else:
                accession = record.get("accession")
                record_id = record.get("id")
                identifier = accession if isinstance(accession, str) and accession else record_id

                if not isinstance(identifier, str) or not identifier:
                    raise ValueError(
                        f"Record {record_number} has no usable accession or ID."
                    )

                clean_identifier = identifier.lstrip(">")
                header = f">{clean_identifier}"

            file.write(header + "\n")
            file.write(record["sequence"] + "\n")
