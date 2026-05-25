from bioseq.alignment.smith_waterman import get_best_scores

def refine_hits(query, top_hits):
    
    for seq in top_hits:

        seq["sw_score"] = get_best_scores(query, seq["sequence"])

    return sorted(top_hits, ley = lambda hit: hit["sw_score"], reverse=True)