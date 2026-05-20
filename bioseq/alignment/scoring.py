def check_match(a, b, match_score = 2, mismatch_score = -1):
    if (a == b):
        return match_score
    return mismatch_score
def normalize_gap(gap):
    return -abs(gap)