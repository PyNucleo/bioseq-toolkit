from Bio.Align import substitution_matrices
from functools import lru_cache

@lru_cache(maxsize=None)
def load_matrix(matrix_name):
    if matrix_name is None:
        return None

    try:
        return substitution_matrices.load(matrix_name)
    except FileNotFoundError as error:
        raise ValueError(f"Unknown substitution matrix: {matrix_name!r}") from error
