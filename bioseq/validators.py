validDNABases = ["A", "T", "G", "C"]
def valid_dna(seq):
    seq = seq.upper()
    return (all(i in validDNABases for i in seq))

def validate_k_and_threshold(k, threshold):
    if not isinstance(k, int) or isinstance(k, bool):
        raise TypeError("k must be an integer.")

    if k < 1:
        raise ValueError("k must be a positive integer.")

    if not isinstance(threshold, int) or isinstance(threshold, bool):
        raise TypeError("threshold must be an integer.")

    if threshold < 1:
        raise ValueError("threshold must be a positive integer.")