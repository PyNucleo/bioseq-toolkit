from bioseq.alignment.smith_waterman import get_best_scores

def rank_database_sequences(query, db):
    """
    This typically follows processing a query against a Database using kmer_search within the pipeline, where it performs multiple local pairwise alignments between the query sequence and each Database sequence

    The alignment score obtained by the alignment of query vs database sequence is stored in a dictionary as items to their corresponding database sequence (Which naturally are the keys)

    Finally, this dictionary is sorted in descending order of the local alignment score obtained from local alignment between the query sequence and the database sequence and returned.
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