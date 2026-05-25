def read_fasta_records(file):
    FASTA = open(file, "r")
    records = []
    current_seq = ''
    record = {"header" : None,
              "sequence" : None}

    for line in FASTA:
        line = line.rstrip()
        if line.startswith(">") and current_seq:
            if len(current_seq) % 3 != 0:
                excess = len(current_seq) % 3
                current_seq = current_seq[0 : len(current_seq) - excess]
            record["sequence"] = current_seq
            records.append(record)
            record = {"header" : line,
                      "sequence" : None}
            current_seq = ''
                
        elif line.startswith(">") and not current_seq:
            record["header"] = line
        else: #Is a sequence, not a header
            current_seq += line

    if len(current_seq) % 3 != 0:
        excess = len(current_seq) % 3
        current_seq = current_seq[0 : len(current_seq) - excess]
    record["sequence"] = current_seq
    records.append(record)
        
    return records

def read_fasta_sequences_only(FILE):

    records = read_fasta_records(FILE)

    sequences = [record["sequence"] for record in records]

    return sequences