from database.database_utils import normalize_database
from bioseq.search.kmer_search import kmer_search
from bioseq.search.kmer_index import multi_query_indexed_search
from bioseq.search.similarity_search import rank_by_shared_kmers
from bioseq.search.refinement import refine_hits
from bioseq.fasta_io import read_fasta_records
from database.database_utils import normalize_database

def search(query, database=None, k=3, threshold=1, top_n_hits=10, refinement=False):
    
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
   
   db = normalize_database(database)

   potential_hits = kmer_search(query, db, k, threshold) #Dictionary of sequence : shared_kmer pairs 

   ranked_hits = rank_by_shared_kmers(query, potential_hits) #Ranks based on SW scores; List of dictionaries, each containing info about one sequence 

   ranked_hits = ranked_hits[:top_n_hits]

   if refinement:
     return refine_hits(query, ranked_hits)
   
   if not refinement:
        return ranked_hits

def search_many(query_fasta, database=None, k=3, threshold=1, top_n_hits=10, indexed=True, refinement=False):
   
    query_records = read_fasta_records(query_fasta)
    
    if indexed:
        db = normalize_database(database)

        hits = multi_query_indexed_search(query_records, db, k, threshold)

        #Sort the hits per each query
        for query_result in hits:
            query_result["indexed_hits"].sort(
                key=lambda hit: (-hit["shared_kmers"], hit["id"])
            )
            query_result["indexed_hits"] = query_result["indexed_hits"][:top_n_hits]


        return hits
    
    else:
        results = []
        for query in query_records:
            results.append(kmer_search(query["sequence"],
                                        database, 
                                        k, 
                                        threshold, 
                                        top_n_hits, 
                                        refinement))

        return results




        


