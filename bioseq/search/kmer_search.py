def kmer_search(query, DB, k, threshold):

    """
    Find candidate database sequences using shared k-mers.

    The query sequence is split into unique words of length k. Each database
    sequence is also split into k-mers, then compared against the query words.
    Sequences with fewer shared k-mers than the threshold are discarded.

    After threshold filtering, candidates are also filtered by a relative score
    cutoff, keeping only sequences whose shared k-mer count is at least a fixed
    fraction of the best candidate's count.

    Parameters
    ----------
    query : str
        Query sequence.
    DB : SequenceDatabase
        Database object containing biological sequences.
    k : int
        K-mer size.
    threshold : int
        Minimum number of shared k-mers required to keep a sequence.

    Returns
    -------
    dict[str, int]
        Dictionary mapping candidate sequences to their shared k-mer counts.
    """
    max_shared_kmers = 0

    DB_SEQUENCES = DB.get_sequences()
    sequences_and_shared_kmers = {x : 0 for x in DB_SEQUENCES}

    query_words = generate_kmers(query, k)

    for db_seq in DB_SEQUENCES:

        temp_shared_kmers = get_shared_kmers(query_words, db_seq, k)

        if temp_shared_kmers < threshold:
            sequences_and_shared_kmers.pop(db_seq)
        
        else:
            if temp_shared_kmers > max_shared_kmers:
                max_shared_kmers = temp_shared_kmers

            sequences_and_shared_kmers[db_seq] = temp_shared_kmers

    return filter_by_relative_score(max_shared_kmers, sequences_and_shared_kmers)


def generate_kmers(seq, k):
    """
    Generate unique k-mers from a sequence.

    Parameters
    ----------
    seq : str
        Input sequence.
    k : int
        K-mer size.

    Returns
    -------
    set[str]
        Set of unique k-mers from the sequence
    """
    return set(seq[i : i + k] for i in range(len(seq) - k + 1))

def get_shared_kmers(query_words, seq, k):
    """
    Count how many query k-mers are present in a database sequence.

    Parameters
    ----------
    query_words : set[str]
        Unique k-mers generated from the query sequence.
    seq : str
        Database sequence to compare against the query words.
    k : int
        K-mer size.

    Returns
    -------
    int
        Number of query k-mers found in the database sequence.
    """
    
    db_seq_words = generate_kmers(seq, k)

    shared_kmers = 0

    for word in query_words:
        if word in db_seq_words:
            shared_kmers += 1

    return shared_kmers

def filter_by_relative_score(max_kmers, candidates_dict, ratio = 0.3):
    """
    Filter candidate sequences by relative shared k-mer count.

    A sequence is kept if its shared k-mer count is at least:
    max_kmers * ratio

    Parameters
    ----------
    max_kmers : int
        Highest shared k-mer count among all candidate sequences.
    candidates_dict : dict[str, int]
        Dictionary mapping sequences to shared k-mer counts.
    ratio : float, optional
        Relative cutoff compared to the best candidate. Default is 0.3.

    Returns
    -------
    dict[str, int]
        Filtered dictionary of candidate sequences and shared k-mer counts.
    """
    filtered = {
    seq_id: count
    for seq_id, count in candidates_dict.items()
    if count >= max_kmers * ratio
}
    return filtered


    