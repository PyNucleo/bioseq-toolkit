from pathlib import Path
from database.database_utils import normalize_database
from bioseq.search.kmer_search import kmer_search, filter_by_relative_score
from bioseq.search.kmer_index import multi_query_indexed_search
from bioseq.search.similarity_search import rank_by_shared_kmers
from bioseq.search.refinement import refine_hits
from bioseq.fasta_io import read_fasta_records
from database.database_utils import normalize_database
from bioseq.validators import validate_k_and_threshold

def search(query, 
           database, 
           k=3, 
           threshold=1, 
           top_n_hits=10, 
           refinement=False, 
           match_score = 1, 
           mismatch_score = -1,
           gap_penalty = -2, 
           matrix = None
           ):
    
   """
    Run a BLAST-like sequence search pipeline.

    The pipeline uses a fast k-mer prefilter to identify candidate database
    sequences that share short words with the query. By default, candidates are
    ranked by the number of shared k-mers. If refinement is enabled, the top
    k-mer hits are re-ranked using Smith-Waterman local alignment.

    Parameters
    ----------
    query : str
        Query biological sequence.
    database : SequenceDatabase | list[str] | str
        Database input. Can be an existing SequenceDatabase, a list of sequence
        strings, or a FASTA file path.
    k : int, optional
        K-mer size used for word matching. Default is 3.
    threshold : int, optional
        Minimum number of shared k-mers required for a database sequence to be
        kept as a candidate. Default is 1.
    top_n_hits : int, optional
        Maximum number of candidate hits to keep before optional refinement.
        Default is 10.
    refinement : bool, optional
        If True, run Smith-Waterman scoring on the top k-mer hits and re-rank
        them by local alignment score. Default is False.

    Returns
    -------
    list[dict]
        Search hits sorted by shared k-mers if refinement is False, or by
        Smith-Waterman score if refinement is True.
    """

   validate_k_and_threshold(k, threshold)

   db = normalize_database(database)

   potential_hits = kmer_search(query, db, k, threshold) #Dictionary of sequence : shared_kmer pairs 

   ranked_hits = rank_by_shared_kmers(query, potential_hits) #Ranks based on SW scores; List of dictionaries, each containing info about one sequence 

   ranked_hits = ranked_hits[:top_n_hits]

   if refinement:
     return refine_hits(query, ranked_hits, gap_penalty, matrix, match_score, mismatch_score)
   
   if not refinement:
        return ranked_hits

def run_indexed_multi_search(query_records, database, k, threshold, top_n_hits, refinement):
    if refinement:
        raise NotImplementedError("Indexed multi-search refinement is not implemented yet.")
    
    db = normalize_database(database)

    results = multi_query_indexed_search(query_records, db, k, threshold)

    for query_result in results:
        hits = query_result["query_hits"]

        if hits:
            max_kmers = max(hit["shared_kmers"] for hit in hits)
            hits = filter_by_relative_score(max_kmers, hits)

        ranked_hits = rank_by_shared_kmers(
            query_result["query_sequence"],
            hits
        )

        query_result["query_hits"] = ranked_hits[:top_n_hits]

    return results

def run_regular_multi_search(query_records, database, k, threshold, top_n_hits, refinement):
    results = []
    for query in query_records:
        temp_result = search(query["sequence"],
                                    database=database, 
                                    k=k, 
                                    threshold=threshold,
                                    top_n_hits=top_n_hits,
                                    refinement=refinement
                                )

        results.append({
            "query_id": query["id"],
            "query_sequence": query["sequence"],
            "query_hits": temp_result
            }
        )

    

    return results

def multi_search(query_fasta,
                 database, 
                 k=3, 
                 threshold=1, 
                 top_n_hits=10, 
                 indexed=True, 
                 refinement=False
                ):

    validate_k_and_threshold(k, threshold)

    if isinstance(query_fasta, (str, Path)):
        query_records = read_fasta_records(query_fasta)
    else:
        query_records = query_fasta
    database = normalize_database(database)

    if indexed:
        return run_indexed_multi_search(query_records, database, k, threshold, top_n_hits, refinement)

    return run_regular_multi_search(query_records, database, k, threshold, top_n_hits, refinement)




        


