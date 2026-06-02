from bioseq.alignment.smith_waterman import get_best_scores
from database.database_utils import normalize_database

def rank_by_shared_kmers(query, db):
    """
    Convert k-mer candidate hits into ranked hit dictionaries.

    This function receives the list of hit dictionaries returned by
    kmer_search() and sorts them by shared k-mer count in descending order.

    Note
    ----
    Despite the function name, this does not currently rank by alignment score.
    Smith-Waterman scoring is handled separately by refine_hits().

    Parameters
    ----------
    query : str
        Query sequence. Currently unused, but kept for pipeline compatibility.
    db : list[dict]
        Candidate hit dictionaries from kmer_search(). Each hit should contain
        "id", "sequence", and "shared_kmers".

    Returns
    -------
    list[dict]
        Hit dictionaries containing:
        - "id": sequence identifier from the normalized database
        - "sequence": database sequence
        - "shared_kmers": number of shared k-mers with the query
    """
    
    ranked_hits = []

    for seq in db:
        
        ranked_hits.append({
            "id" : seq["id"],
            "sequence": seq["sequence"],
            "shared_kmers": seq["shared_kmers"]
        })
    
    ranked_hits = sorted(ranked_hits, key = lambda hit: hit["shared_kmers"], reverse=True)

    return ranked_hits
