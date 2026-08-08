validDNABases = ["A", "T", "G", "C"]


def report_invalid_symbols_and_positions(seq):
    seq = seq.upper()

    invalid_positions = []
    invalid_symbols = set()

    for position, symbol in enumerate(seq, start=1):
        if symbol not in validDNABases:
            invalid_positions.append(position)
            invalid_symbols.add(symbol)


    report = {
        "invalid_positions": invalid_positions,
        "invalid_symbols": sorted(invalid_symbols)
    }

    return report

def valid_dna(seq):
    seq = seq.upper()
    return all(i in validDNABases for i in seq)

def validate_k_and_threshold(k, threshold):
    if not isinstance(k, int) or isinstance(k, bool):
        raise TypeError("k must be an integer.")

    if k < 1:
        raise ValueError("k must be a positive integer.")

    if not isinstance(threshold, int) or isinstance(threshold, bool):
        raise TypeError("threshold must be an integer.")

    if threshold < 1:
        raise ValueError("threshold must be a positive integer.")