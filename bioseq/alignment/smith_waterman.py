from .scoring import score_pair, normalize_gap, build_scoring_metadata
from .alignment_stats import get_alignment_stats

def initialize_grid(s1, s2):
    grid = [[0] * (len(s2) + 1) for _ in range(len(s1) + 1)]
    return grid
def initialize_movementGrid(s1, s2):
    tempMoveGrid = [[[None]] * (len(s2) + 1) for _ in range(len(s1) + 1)]
    return tempMoveGrid
def dynamic_movement_filling(d, l, v):
    
    tempArray = []
    temp = max(d, l, v, 0)
    
    if (temp == d):
        tempArray.append("diagonal")
        
    if (temp == l):
        tempArray.append("left")
    
    if (temp == v):
        tempArray.append("vertical")
    
    if (d < 0 and l < 0 and v < 0):
        tempArray.append(None)
    
    return tempArray

def fill_local(s1, s2, matrix, grid,movement_grid, gap_penalty, match_s=2, mismatch_s=-1):
   
    gap_penalty = normalize_gap(gap_penalty)

    
    tempMax_At = []
    tempMax = 0
    
    for i in range(1, len(s1) + 1):
        
        for j in range(1, len(s2) + 1):
            match_score = score_pair(s1[i - 1], s2[j - 1], matrix, match_score=match_s, mismatch_score=mismatch_s)
            
            diagonal = grid[i - 1][j - 1] + match_score
            horizontal = grid[i][j - 1] + gap_penalty
            vertical = grid[i - 1][j] + gap_penalty 

            
            score = max(diagonal, horizontal, vertical, 0)
            
            if (score > tempMax):
                tempMax = score
                tempMax_At = [(i, j)]
            elif (score == tempMax and score > 0):
                tempMax_At.append((i, j))
                
            movement_grid[i][j] = dynamic_movement_filling(diagonal, horizontal, vertical)
            grid[i][j] = score
    return tempMax, tempMax_At

def local_trace(movementMatrix, s1, s2, row, col, algn1, algn2, return_all, scoring_matrix):
    if (scoring_matrix[row][col] == 0):
        algn1, algn2 = algn1[::-1], algn2[::-1]
        return [(algn1, algn2)]
    
    result = []
    
    if return_all:
        
        for move in movementMatrix[row][col]:
            
            if (move == "diagonal"):
                new_branch1 = algn1 + s1[row - 1]
                new_branch2 = algn2 + s2[col - 1]
                result.extend(local_trace(movementMatrix, s1, s2, row - 1, col - 1, new_branch1, new_branch2, return_all, scoring_matrix))
            elif (move == "left"):
                new_branch1 = algn1 + "-"
                new_branch2 = algn2 + s2[col - 1]
                result.extend(local_trace(movementMatrix, s1, s2, row, col - 1, new_branch1, new_branch2, return_all, scoring_matrix))
            elif (move == "vertical"):
                new_branch1 = algn1 + s1[row - 1]
                new_branch2 = algn2 + "-"
                result.extend(local_trace(movementMatrix, s1, s2, row - 1, col, new_branch1, new_branch2, return_all, scoring_matrix))
    else:
            move = movementMatrix[row][col][0]
            
            if (move == "diagonal"):
                new_branch1 = algn1 + s1[row - 1]
                new_branch2 = algn2 + s2[col - 1]
                branch_align = local_trace(movementMatrix, s1, s2, row - 1, col - 1, new_branch1, new_branch2, return_all, scoring_matrix)
            elif (move == "left"):
                new_branch1 = algn1 + "-"
                new_branch2 = algn2 + s2[col - 1]
                branch_align = local_trace(movementMatrix, s1, s2, row, col - 1, new_branch1, new_branch2, return_all, scoring_matrix)
            elif (move == "vertical"):
                new_branch1 = algn1 + s1[row - 1]
                new_branch2 = algn2 + "-"
                branch_align = local_trace(movementMatrix, s1, s2, row - 1, col, new_branch1, new_branch2, return_all, scoring_matrix)
            result.extend(branch_align)
    
    return result

def get_best_scores(s1, s2, gap_penalty, matrix=None, match=2, mismatch=-1):

    gap_penalty = normalize_gap(gap_penalty)
    
    grid = initialize_grid(s1, s2)
    
    tempMax_At = []
    tempMax = 0
    
    for i in range(1, len(s1) + 1):
        
        for j in range(1, len(s2) + 1):
            match_score = score_pair(s1[i - 1], s2[j - 1], matrix, match_score=match, mismatch_score=mismatch)
            
            diagonal = grid[i - 1][j - 1] + match_score
            horizontal = grid[i][j - 1] + gap_penalty
            vertical = grid[i - 1][j] + gap_penalty 

            
            score = max(diagonal, horizontal, vertical, 0)
            
            if (score > tempMax):
                tempMax = score
                tempMax_At = [(i, j)]
            elif (score == tempMax and score > 0):
                tempMax_At.append((i, j))
            grid[i][j] = score
            
    return tempMax, tempMax_At
def local_alignment(s1, s2, match=2, mismatch=-1, gap_penalty=-2, matrix=None, return_all=False, structured=True):
    if (s1 == "" or s2 == ""):
        raise ValueError("Sequence(s) cannot be empty.")
    
    gap_penalty = normalize_gap(gap_penalty)

    grid = initialize_grid(s1, s2)
    
    movement_grid = initialize_movementGrid(s1, s2)
    
    best_score, start_at = fill_local(s1, s2, matrix, grid,movement_grid, gap_penalty, match_s=match,mismatch_s=mismatch)

    all_alignments = []
    
    if not return_all:
        start_at = start_at[:1]

    for pos in start_at:
        all_alignments.extend(local_trace(movement_grid, s1, s2, pos[0], pos[1], "", "", return_all, grid))

    
    if structured:
        temp_structure = {
            "algorithm": "Smith-Waterman",
            "mode": "local",
            "sequence_1": s1,
            "sequence_2": s2,
            "score": best_score,
            "scoring":{
                **build_scoring_metadata(match, mismatch, gap_penalty, matrix) 
            },
            "best_positions": start_at,
            "num_alignments": len(all_alignments),
            "alignments":[]
        }
        
        for i in range(len(all_alignments)):
            algn_1 = all_alignments[i][0]
            algn_2 = all_alignments[i][1]

            alignment_stats = get_alignment_stats(algn_1, algn_2)

            temp_structure["alignments"].append({
                "aligned_sequence_1": algn_1,
                "aligned_sequence_2": algn_2,
                **alignment_stats
            })
        
        return temp_structure
    
    
    
    
    return all_alignments
    





