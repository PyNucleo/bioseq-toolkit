validDNABases = ["A", "T", "G", "C"]
def validate_dna(seq):
    seq = seq.upper()
    return (all(i in validDNABases for i in seq))
