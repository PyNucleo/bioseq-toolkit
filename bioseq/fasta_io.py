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

            record = {
                "header": line,
                "sequence": ""
            }

            current_seq = ""

        else:

            current_seq += line

    if record is not None:

        record["sequence"] = current_seq
        records.append(record)

    return records

def read_fasta_sequences_only(FILE):

    records = read_fasta_records(FILE)

    sequences = [record["sequence"] for record in records]

    return sequences