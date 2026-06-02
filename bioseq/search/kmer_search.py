def kmer_search(query, db, k, threshold):

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
    is_valid_query_length = validate_kmer_params(query, k, threshold)

    if not is_valid_query_length:
        return []
    
    max_shared_kmers = 0

    db_records = db.get_sequences() #List of Dicts

    results = []

    query_words = generate_kmers(query, k)

    for db_seq in db_records:

        temp_shared_kmers = get_shared_kmers(query_words, db_seq["sequence"], k)

        if temp_shared_kmers >= threshold:
            results.append({"id": db_seq["id"],
                            "sequence": db_seq["sequence"],
                            "shared_kmers": temp_shared_kmers})
        
        if temp_shared_kmers > max_shared_kmers:
            max_shared_kmers = temp_shared_kmers


    return filter_by_relative_score(max_shared_kmers, results)

def validate_kmer_params(seq, k, threshold=None):
    if not isinstance(k, int):
        raise TypeError("k must be an integer.")

    if k <= 0:
        raise ValueError("k must be a positive integer.")

    if k > len(seq):
        return False

    if threshold is not None:
        if not isinstance(threshold, int):
            raise TypeError("threshold must be an integer.")

        if threshold < 0:
            raise ValueError("threshold must be greater than or equal to 0.")

    return True

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

    is_valid_length = validate_kmer_params(seq, k)

    if not is_valid_length:
        return set()

    seq = seq.upper()

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
    filtered = [
    d for d in candidates_dict
    if d["shared_kmers"] >= max_kmers * ratio
]
    return filtered


    