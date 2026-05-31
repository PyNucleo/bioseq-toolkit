from .scoring import check_match, normalize_gap
from .alignment_stats import get_identity, get_matches_mismatches_gaps, get_alignment_stats

def initialize_grid(s1, s2, gap_penalty):
    grid = [[0] * (len(s2) + 1) for _ in range(len(s1) + 1)]
    rows = len(s1) + 1
    cols = len(s2) + 1
    
    for i in range(1, rows):
        grid[i][0] = grid[i - 1][0] + gap_penalty
    
    for j in range(1, cols):
        grid[0][j] = grid[0][j - 1] + gap_penalty
        
    return grid

def initialGlobalMovementMatrix(s1, s2):
    movements = [[None] * (len(s2) + 1) for _ in range(len(s1) + 1)]
    
    
    for col in range(1, len(s2) + 1):
        movements[0][col] = ["left"]
    for row in range(1, len(s1) + 1):
        movements[row][0] = ["up"]
    return movements
def fill_matrix(grid, s1, s2, match_score, mismatch_score, gap_penalty):
    for row in range(1, len(grid)):
        for column in range(1, len(s2) + 1):
            a, b = s1[row - 1], s2[column - 1]
            match_case = check_match(a, b, match_score, mismatch_score)
            grid[row][column] = max((grid[row - 1][column - 1] + match_case), grid[row][column - 1] + gap_penalty, grid[row - 1][column] + gap_penalty)
    return grid
        
def dynamicMovementMatrix(grid, row, column, check_match, gap):

    movement_case = []
    diagonal = grid[row - 1][column - 1] + check_match
    vertical_up = grid[row - 1][column] + gap
    horizontal_left = grid[row][column - 1] + gap
    
    max_score = max(diagonal, vertical_up, horizontal_left)
    
    if (diagonal == max_score):
        movement_case.append("diag")
        
    if (horizontal_left == max_score):
        movement_case.append("left")
    
    if (vertical_up == max_score):
        movement_case.append("up")
    
    return movement_case
def compute_movement_matrix(s1, s2, grid, match_score, mismatch_score, gap):
    movements = initialGlobalMovementMatrix(s1, s2)
    for row in range(1, len(s1) + 1):
        for column in range(1, len(s2) + 1):
            a, b = s1[row - 1], s2[column - 1]
            match_case = check_match(a, b, match_score, mismatch_score)
            movements[row][column] = dynamicMovementMatrix(grid, row, column, match_case, gap)
    return movements

def trace(movements, s1, s2, row, col, algn1, algn2, return_all):
    if row == 0 and col == 0:
        algn1, algn2 = algn1[::-1], algn2[::-1]
        return [(algn1, algn2)]
    
    result = []
    
    if return_all:
        for move in movements[row][col]:
            
            if move == "diag":
                new_branch1 = algn1 + s1[row - 1]
                new_branch2 = algn2 + s2[col - 1]
                branch_align = trace(movements, s1, s2, row - 1, col - 1, new_branch1, new_branch2, return_all)
            elif move == "up":
                new_branch1 = algn1 + s1[row - 1]
                new_branch2 = algn2 + "-"
                branch_align = trace(movements, s1, s2, row - 1, col, new_branch1, new_branch2, return_all)
            else:
                new_branch1 = algn1 + "-"
                new_branch2 = algn2 + s2[col - 1]
                branch_align = trace(movements, s1, s2, row, col - 1, new_branch1, new_branch2, return_all)
    
    else:
            move = movements[row][col][0]
            
            if move == "diag":
                new_branch1 = algn1 + s1[row - 1]
                new_branch2 = algn2 + s2[col - 1]
                branch_align = trace(movements, s1, s2, row - 1, col - 1, new_branch1, new_branch2, return_all)
            elif move == "up":
                new_branch1 = algn1 + s1[row - 1]
                new_branch2 = algn2 + "-"
                branch_align = trace(movements, s1, s2, row - 1, col, new_branch1, new_branch2, return_all)
            else:
                new_branch1 = algn1 + "-"
                new_branch2 = algn2 + s2[col - 1]
                branch_align = trace(movements, s1, s2, row, col - 1, new_branch1, new_branch2, return_all)

    result.extend(branch_align)
    
    return result

def global_alignment(s1, s2, match=1, mismatch=-1, gap_penalty=-2, return_all = False, structured=True):
    gap_penalty = normalize_gap(gap_penalty)
        
    # 1. Initialize grid
    grid = initialize_grid(s1, s2, gap_penalty)

    # 2. Fill scoring matrix
    grid = fill_matrix(grid, s1, s2, match, mismatch, gap_penalty)

    # 3. Compute movement matrix
    movements = compute_movement_matrix(s1, s2, grid, match, mismatch, gap_penalty)
        
    # 4. Traceback (ALL alignments)
    alignments = trace(movements, s1, s2, len(s1), len(s2), "", "", return_all)

    if structured:
        temp_structure = {
            "algorithm": "Needleman-Wunsch",
            "mode": "global",
            "sequence_1": s1,
            "sequence_2": s2,
            "score": grid[len(s1)][len(s2)],
            "scoring":{
                "match": match,
                "mismatch": mismatch,
                "gap": gap_penalty,
                "matrix": None,
                "gap_model": "linear"  
            },
            "num_alignments": len(alignments),
            "alignments":[]
        }
        
        for i in range(len(alignments)):
            algn_1 = alignments[i][0]
            algn_2 = alignments[i][1]

            matches, mismatches, gaps, gap_columns = get_matches_mismatches_gaps(algn_1, algn_2)

            alignment_stats = get_alignment_stats(algn_1, algn_2)

            temp_structure["alignments"].append({
                "aligned_sequence_1": algn_1,
                "aligned_sequence_2": algn_2,
                **alignment_stats
            })


        return temp_structure
    
    else:
        return alignments
