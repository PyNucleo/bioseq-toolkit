

def get_identity(s1, s2, including_gaps = True):
    if len(s1) != len(s2):
        raise ValueError("Aligned sequences must have the same length.")
    
    tot = 0
    l = len(s1)

    for position in range(l):
            
            if (s1[position] == s2[position]):
                tot += 1
    
    if including_gaps:
        return round(tot / l, 3)
    
    length_without_gaps = l - (s1.count("-") + s2.count("-")) 

    if length_without_gaps == 0:
         return 0
    
    return round(tot / length_without_gaps, 3)

def get_matches(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Aligned sequences must have the same length.")
    
    l = len(s1)
    matches = 0

    for position in range(l):
            if s1[position] == s2[position] and s1[position] != "-":
                matches += 1
    return matches

def get_mismatches(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Aligned sequences must have the same length.")
    
    l = len(s1)
    mismatches = 0

    for position in range(l):
            if (s1[position] != s2[position] and s1[position] != "-" and s2[position] != "-"):
                mismatches += 1
    return mismatches

def get_gaps(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Aligned sequences must have the same length.")
    
    l = len(s1)
    gaps = 0

    for position in range(l):
            if (s1[position] == "-"):
                gaps += 1
            
            if (s2[position] == "-"):
                gaps += 1
    
    return gaps

def get_gap_columns(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Aligned sequences must have the same length.")
    
    gap_columns = 0

    for i in range(len(s1)):
         if s1[i] == "-":
              gap_columns += 1
              continue
         elif s2[i] == "-":
              gap_columns += 1
              continue
    
    return gap_columns

def get_matches_mismatches_gaps(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Aligned sequences must have the same length.")
    
    l = len(s1)
    matches = mismatches = gaps = gap_columns = 0

    for position in range(l):
            if (s1[position] == s2[position] and (s1[position] != "-") and (s2[position] != "-")):
                matches += 1
            
            else:

                if (s1[position] != "-") and (s2[position] != "-"):
                     mismatches += 1
                
                else:
                    
                    if (s1[position] == "-") and (s2[position] == "-"):
                         gaps += 2
                    else:
                         gaps += 1
                    gap_columns += 1
    
    return matches, mismatches, gaps, gap_columns

def get_alignment_stats(algn_1, algn_2):
     
    matches, mismatches, gaps, gap_columns = get_matches_mismatches_gaps(algn_1, algn_2)

    temp_structure = {
    "alignment_length": len(algn_1),
    "matches": matches,
    "mismatches": mismatches,
    "gaps": gaps,
    "gap_columns": gap_columns,
    "identity": get_identity(algn_1, algn_2),
    "identity_excluding_gaps": get_identity(algn_1, algn_2, including_gaps=False),
    "similarity": None
    }

    return temp_structure