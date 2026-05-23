from bioseq.alignment.smith_waterman import get_best_scores

def rank_database_sequences(query, DB):
    
    scoresDict = dict()

    DB_SEQUENCES = DB.get_sequences()

    for seq in DB_SEQUENCES:

        score, _ = get_best_scores(query, seq, 2)
        
        
        scoresDict[seq] = score
    
    scoresDict = sorted(scoresDict.items(), key = lambda x: x[1], reverse=True)

    return(scoresDict)