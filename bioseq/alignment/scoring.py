from .substitution_matrices import load_matrix

def get_scoring_matrix(matrix=None):
    if matrix is None:
        return None

    if isinstance(matrix, str):
        return load_matrix(matrix)

    return matrix

def score_pair(a, b, loaded_matrix = None, match_score = 2, mismatch_score = -1):
    a = a.upper()
    b = b.upper()

    if loaded_matrix is None:
        return match_score if a == b else mismatch_score

    loaded_matrix = get_scoring_matrix(loaded_matrix)

    if a not in loaded_matrix.alphabet or b not in loaded_matrix.alphabet:
        raise ValueError(f"Residue pair ({a!r}, {b!r}) is not supported by the selected substitution matrix.") 
    return loaded_matrix[a, b]


def normalize_gap(gap):
    return -abs(gap)

def build_scoring_metadata(match, mismatch, gap_penalty, matrix=None):
    if matrix == None:
        return {
            "match": match,
            "mismatch": mismatch,
            "gap_penalty": gap_penalty,
            "matrix": None,
            "gap_model": "linear"
        }
    else:
        return {
            "gap_penalty": gap_penalty,
            "matrix": matrix,
            "gap_model": "linear"
        }
