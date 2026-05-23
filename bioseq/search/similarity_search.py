from bioseq.alignment.smith_waterman import get_best_scores

def rank_database_sequences(query, DB):
    """
    This typically follows processing a query against a Database using kmer_search within the pipeline, where it performs multiple local pairwise alignments between the query sequence and each Database sequence

    The alignment score obtained by the alignment of query vs database sequence is stored in a dictionary as items to their corresponding database sequence (Which naturally are the keys)

    Finally, this dictionary is sorted in descending order of the local alignment score obtained from local alignment between the query sequence and the database sequence and returned.
    """
    scores_dict = dict()

    DB_SEQUENCES = DB.get_sequences()

    for seq in DB_SEQUENCES:

        score, _ = get_best_scores(query, seq, 2)
        
        
        scores_dict[seq] = score
    
    scores_dict = sorted(scores_dict.items(), key = lambda x: x[1], reverse=True)

    return(scores_dict)