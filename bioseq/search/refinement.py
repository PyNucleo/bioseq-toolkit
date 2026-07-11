from bioseq.alignment.smith_waterman import get_best_scores


def refine_hits(query, top_hits, gap_penalty=-2, matrix=None, match_score=1, mismatch_score=-1):
    """
    Refine k-mer search hits using Smith-Waterman local alignment.

    The k-mer search step is fast but approximate: it ranks candidates by
    shared words rather than by an actual alignment score. This function takes
    those candidate hits, computes a Smith-Waterman score for each one, stores
    the score and best-scoring matrix positions, then returns the hits sorted
    by local alignment score.

    Parameters
    ----------
    query : str
        The query sequence.
    top_hits : list[dict]
        Candidate hits from the k-mer search pipeline. Each hit is expected to
        contain at least a "sequence" key.
    gap_penalty : int, optional
        Gap penalty used by Smith-Waterman. Positive values are converted to
        negative values by the scoring utility. Default is -2.
    matrix : str, optional
        The matrix to be used for scoring. If not passed, linear match/mismatch
        scoring will be used. If an invalid matrix is passed, raises "KeyError".
    match_score : int, optional
        The score used for matching residues if no matrix is passed.
    mismatch_score : int, optional
        The score used for mismatching residues if no matrix is passed.

    Returns
    -------
    list[dict]
        Refined hit dictionaries sorted by descending Smith-Waterman score.
        Each returned hit includes:
        - "sw_score": best local alignment score
        - "best_positions": matrix positions where the best score occurs
    """
    refined_hits = []

    for hit in top_hits:
        refined_hit = hit.copy()

        score, positions = get_best_scores(
            query,
            refined_hit["sequence"],
            gap_penalty,
            matrix=matrix,
            match=match_score,
            mismatch=mismatch_score
        )

        refined_hit["sw_score"] = score
        refined_hit["best_positions"] = positions

        refined_hits.append(refined_hit)

    return sorted(
        refined_hits,
        key=lambda hit: hit["sw_score"],
        reverse=True
    )