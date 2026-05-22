from .scoring import check_match, normalize_gap

def initialize_grid(s1, s2):
    grid = [[0] * (len(s2) + 1) for _ in range(len(s1) + 1)]
    return grid
def initialize_movementGrid(s1, s2):
    tempMoveGrid = [[[None]] * (len(s2) + 1) for _ in range(len(s1) + 1)]
    return tempMoveGrid
def dynamic_movement_filling(d, l, v, ):
    
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

def fill_local(s1, s2, grid,movement_grid, gap_penalty):
   
    gap_penalty = normalize_gap(gap_penalty)

    
    tempMax_At = []
    tempMax = 0
    
    for i in range(1, len(s1) + 1):
        
        for j in range(1, len(s2) + 1):
            match_score = check_match(s1[i - 1], s2[j - 1])
            
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
                branch_align = local_trace(movementMatrix, s1, s2, row - 1, col - 1, new_branch1, new_branch2, return_all, scoring_matrix)
            elif (move == "left"):
                new_branch1 = algn1 + "-"
                new_branch2 = algn2 + s2[col - 1]
                branch_align = local_trace(movementMatrix, s1, s2, row, col - 1, new_branch1, new_branch2, return_all, scoring_matrix)
            elif (move == "vertical"):
                new_branch1 = algn1 + s1[row - 1]
                new_branch2 = algn2 + "-"
                branch_align = local_trace(movementMatrix, s1, s2, row - 1, col, new_branch1, new_branch2, return_all, scoring_matrix)
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

def get_best_scores(s1, s2, gap_penalty):
    
    gap_penalty = normalize_gap(gap_penalty)
    
    grid = initialize_grid(s1, s2)
    
    tempMax_At = []
    tempMax = 0
    
    for i in range(1, len(s1) + 1):
        
        for j in range(1, len(s2) + 1):
            match_score = check_match(s1[i - 1], s2[j - 1])
            
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
def local_alignment(s1, s2, return_all = False):
    
    grid = initialize_grid(s1, s2)
    
    movement_grid = initialize_movementGrid(s1, s2)
    
    best_Score, startAt = fill_local(s1, s2, grid,movement_grid, -2)

    all_alignments = []
    
    for pos in startAt:
        all_alignments.append(local_trace(movement_grid, s1, s2, pos[0], pos[1], "", "", return_all, grid))
    
    return all_alignments
    





