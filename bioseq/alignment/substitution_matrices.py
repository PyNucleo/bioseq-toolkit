from Bio.Align import substitution_matrices
from functools import lru_cache

@lru_cache(maxsize=None)
def load_matrix(matrix_name):
    if matrix_name is None:
        return None
    
    return substitution_matrices.load(matrix_name)
