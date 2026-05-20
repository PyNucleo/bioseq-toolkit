from bioseq.alignment.smith_waterman import get_best_scores

def rank_database_sequences(query):
    
    database = {
    "seq1": "GATTACA",
    "seq2": "GATTTCA",
    "seq3": "TTACAGG",
    "seq4": "CCCCCCC",
    "seq5": "GACTATA",
    "seq6": "GATACCA",
    "seq7": "TTACTAA",
    "seq8": "GGGGATTACAAGG",
    "seq9": "ATCGATCG",
    "seq10": "GATTGCA"
}
    
    scoresDict = dict()
    for seq in database:

        score, _ = get_best_scores(query, database[seq], 2)
        
        
        scoresDict[seq] = score
    
    scoresDict = sorted(scoresDict.items(), key = lambda x: x[1], reverse=True)

    return(scoresDict)