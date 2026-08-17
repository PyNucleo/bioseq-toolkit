from pathlib import Path
from database.database_utils import normalize_database
from bioseq.search.kmer_search import kmer_search, filter_by_relative_score
from bioseq.search.kmer_index import multi_query_indexed_search
from bioseq.search.similarity_search import rank_by_shared_kmers
from bioseq.search.refinement import refine_hits
from bioseq.fasta_io import read_fasta_records
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
    match_score : int, optional
        Match score used during refinement when ``matrix`` is None. Ignored
        when a substitution matrix is selected. Default is 1.
    mismatch_score : int, optional
        Mismatch score used during refinement when ``matrix`` is None. Ignored
        when a substitution matrix is selected. Default is -1.
    gap_penalty : int, optional
        Linear gap penalty used during refinement in simple and matrix modes.
        Default is -2.
    matrix : str, optional
        Optional substitution-matrix name. When supplied, matrix residue-pair
        scores replace simple match/mismatch scoring. Default is None.

    All scoring parameters are inactive when ``refinement`` is False.

    Returns
    -------
    list[dict]
        Search hits sorted by shared k-mers if refinement is False, or by
        Smith-Waterman score if refinement is True.
    """

   validate_k_and_threshold(k, threshold)

   db = normalize_database(database)

   potential_hits = kmer_search(query, db, k, threshold) # Candidate hit dictionaries after k-mer filtering.

   # Rank k-mer candidates by descending shared-k-mer count, then ascending ID.
   ranked_hits = rank_by_shared_kmers(query, potential_hits)

   ranked_hits = ranked_hits[:top_n_hits]

   if refinement:
     return refine_hits(
         query,
         ranked_hits,
         gap_penalty=gap_penalty,
         matrix=matrix,
         match_score=match_score,
         mismatch_score=mismatch_score,
     )
   
   if not refinement:
        return ranked_hits

def run_indexed_multi_search(
    query_records,
    database,
    k,
    threshold,
    top_n_hits,
    refinement,
    match_score=1,
    mismatch_score=-1,
    gap_penalty=-2,
    matrix=None,
):
    """Run indexed candidate selection and optionally refine each query's hits.

    Refinement uses the same ``refine_hits()`` contract as regular search and
    runs only after the existing candidate filtering, ranking, and
    ``top_n_hits`` selection. Scoring parameters are inactive when refinement
    is disabled. In matrix mode, ``match_score`` and ``mismatch_score`` are
    accepted but ignored; ``gap_penalty`` remains active.
    """
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

        selected_hits = ranked_hits[:top_n_hits]

        if refinement:
            selected_hits = refine_hits(
                query_result["query_sequence"],
                selected_hits,
                gap_penalty=gap_penalty,
                matrix=matrix,
                match_score=match_score,
                mismatch_score=mismatch_score,
            )

        query_result["query_hits"] = selected_hits

    return results

def run_regular_multi_search(
    query_records,
    database,
    k,
    threshold,
    top_n_hits,
    refinement,
    match_score=1,
    mismatch_score=-1,
    gap_penalty=-2,
    matrix=None,
):
    """Run each query through ``search()`` with explicit scoring controls.

    ``search()`` performs optional refinement after its established candidate
    selection and ``top_n_hits`` boundary. Scoring parameters are inactive
    when refinement is disabled. In matrix mode, simple match/mismatch values
    are accepted but ignored, while the linear gap penalty remains active.
    """
    results = []
    for query in query_records:
        temp_result = search(
            query["sequence"],
            database=database,
            k=k,
            threshold=threshold,
            top_n_hits=top_n_hits,
            refinement=refinement,
            match_score=match_score,
            mismatch_score=mismatch_score,
            gap_penalty=gap_penalty,
            matrix=matrix,
        )

        results.append({
            "query_id": query["id"],
            "query_sequence": query["sequence"],
            "query_hits": temp_result
            }
        )

    return results

def multi_search(
    query_fasta,
    database,
    k=3,
    threshold=1,
    top_n_hits=10,
    indexed=True,
    refinement=False,
    match_score=1,
    mismatch_score=-1,
    gap_penalty=-2,
    matrix=None,
):
    """Search multiple queries using regular or indexed candidate discovery.

    Parameters
    ----------
    query_fasta : str | pathlib.Path | list[dict]
        Query FASTA path or query records containing ``id`` and ``sequence``.
    database : SequenceDatabase | list[str] | str
        Database input accepted by ``normalize_database()``.
    k : int, optional
        K-mer size. Default is 3.
    threshold : int, optional
        Minimum shared-k-mer count. Default is 1.
    top_n_hits : int, optional
        Maximum selected candidates per query before optional refinement.
        Default is 10.
    indexed : bool, optional
        Use the indexed multi-query candidate path when True and the regular
        ``search()`` path when False. Default is True.
    refinement : bool, optional
        Refine each query's selected candidates with Smith-Waterman. Supported
        in both regular and indexed modes. Default is False.
    match_score : int, optional
        Match score used during refinement when ``matrix`` is None. Ignored
        when a substitution matrix is selected. Default is 1.
    mismatch_score : int, optional
        Mismatch score used during refinement when ``matrix`` is None. Ignored
        when a substitution matrix is selected. Default is -1.
    gap_penalty : int, optional
        Linear gap penalty used during refinement in simple and matrix modes.
        Default is -2.
    matrix : str, optional
        Optional substitution-matrix name. When supplied, matrix residue-pair
        scores replace simple match/mismatch scoring. Default is None.

    Notes
    -----
    All scoring parameters are inactive when ``refinement`` is False. Indexed
    refinement keeps the indexed path's existing selected candidate set and
    applies the same ``refine_hits()`` contract used by regular refinement.
    """

    validate_k_and_threshold(k, threshold)

    if isinstance(query_fasta, (str, Path)):
        query_records = read_fasta_records(query_fasta)
    else:
        query_records = query_fasta
    database = normalize_database(database)

    if indexed:
        return run_indexed_multi_search(
            query_records,
            database,
            k,
            threshold,
            top_n_hits,
            refinement,
            match_score=match_score,
            mismatch_score=mismatch_score,
            gap_penalty=gap_penalty,
            matrix=matrix,
        )

    return run_regular_multi_search(
        query_records,
        database,
        k,
        threshold,
        top_n_hits,
        refinement,
        match_score=match_score,
        mismatch_score=mismatch_score,
        gap_penalty=gap_penalty,
        matrix=matrix,
    )
