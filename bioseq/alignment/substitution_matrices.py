from Bio.Align import substitution_matrices

def load_matrix(matrix_name):
    if matrix_name is None:
        return None
    
    return substitution_matrices.load(matrix_name)
