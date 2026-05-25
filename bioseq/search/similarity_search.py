from bioseq.alignment.smith_waterman import get_best_scores

def rank_by_shared_kmers(query, db):
    """
    Convert k-mer candidate hits into ranked hit dictionaries.

    This function receives the dictionary returned by kmer_search(), where each
    key is a database sequence and each value is the number of k-mers shared
    with the query. It converts those entries into a list of hit dictionaries
    and sorts them by shared k-mer count in descending order.

    Note
    ----
    Despite the function name, this does not currently rank by alignment score.
    Smith-Waterman scoring is handled separately by refine_hits().

    Parameters
    ----------
    query : str
        Query sequence. Currently unused, but kept for pipeline compatibility.
    db : dict[str, int]
        Dictionary mapping database sequences to shared k-mer counts.

    Returns
    -------
    list[dict]
        Hit dictionaries containing:
        - "seq_id": currently None
        - "sequence": database sequence
        - "shared_kmers": number of shared k-mers with the query
    """
    
    ranked_hits = []

    for seq in db:
        
        ranked_hits.append({
            "seq_id" : None,
            "sequence":seq,
            "shared_kmers": db[seq]
        })
    
    ranked_hits = sorted(ranked_hits, key = lambda hit: hit["shared_kmers"], reverse=True)

    return(ranked_hits)