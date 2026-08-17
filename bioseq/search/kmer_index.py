from bioseq.search.kmer_search import generate_kmers


def multi_query_indexed_search(queries, db, k, threshold):
    """
    Search multiple query sequences against a database using a precomputed k-mer index.

    This function is intended for the case where many query sequences are searched
    against the same database. Instead of scanning every database sequence separately
    for every query, the database is indexed once with `index_database_words()`.
    Each query is then broken into k-mers, and those k-mers are looked up directly
    in the database-wide index.

    Parameters
    ----------
    queries : list[dict]
        Normalized query records containing "id" and "sequence".
    db : SequenceDatabase
        Normalized database object containing records with unique IDs.
    k : int
        Previously validated positive k-mer size.
    threshold : int
        Previously validated positive minimum shared-k-mer count.

    Returns
    -------
    list of dict
        One result dictionary per query sequence. Each result has:

        - "query_id": ID of the query sequence
        - "query_sequence": original query sequence
        - "query_hits": list of database hits with shared k-mer counts

        Example:
        [
            {
                "query_id": "query1",
                "query_sequence": "ATGCG",
                "query_hits": [
                    {"id": "seq1", "sequence": "ATGCGT", "shared_kmers": 3},
                    {"id": "seq2", "sequence": "ATGCGA", "shared_kmers": 3}
                ]
            }
        ]

    Notes
    -----
    The database index is built once and reused for all queries. This is the main
    advantage of indexed search when running multiple queries against the same
    database.
    """

    indexed_kmers_structure = index_database_words(db, k)

    query_indexed_search_results = []

    for query in queries:

        query_words = generate_kmers(query["sequence"], k)

        db_sequences_shared_kmer_counts = get_kmer_occurrences(
        query_words,
        indexed_kmers_structure["kmer_index"],
        indexed_kmers_structure["sequence_lookup"],
        threshold
        )

        query_indexed_search_results.append({
            "query_id": query["id"],
            "query_sequence": query["sequence"],
            "query_hits": db_sequences_shared_kmer_counts
        })

    return query_indexed_search_results


def get_kmer_occurrences(words, indexed_kmers, sequence_lookup, threshold):
    """
    Count how many query k-mers are shared with each indexed database sequence.

    For a single query sequence, this function receives the query's k-mers and
    looks each one up in a precomputed database index. Every time a database
    sequence contains one of the query k-mers, that database sequence's
    `shared_kmers` count is incremented.

    Parameters
    ----------
    words : set[str]
        Unique query k-mers.
    indexed_kmers : dict[str, set[str]]
        Inverted k-mer index.
    sequence_lookup : dict[str, str]
        Database ID-to-sequence mapping.
    threshold : int
        Previously validated positive minimum shared-k-mer count.

    Returns
    -------
    list of dict
        Database hit records with shared k-mer counts.

        Example:
        [
            {"id": "seq1", "sequence": "ATGCGT", "shared_kmers": 3},
            {"id": "seq2", "sequence": "ATGAAA", "shared_kmers": 1}
        ]

    Notes
    -----
    The function uses a dictionary internally so each database sequence ID can be
    updated in average constant time. The final result is converted to a list of
    dictionaries because ranked search outputs are easier to sort, filter, and
    display in list form.
    """

    hit_records = {}

    for word in words:

        if word not in indexed_kmers:
            continue

        for db_id in indexed_kmers[word]:

            if db_id not in hit_records:
                hit_records[db_id] = {
                    "id": db_id,
                    "sequence": sequence_lookup[db_id],
                    "shared_kmers": 0
                }

            hit_records[db_id]["shared_kmers"] += 1

    indexed_hits = [
        hit
        for hit in hit_records.values()
        if hit["shared_kmers"] >= threshold
    ]

    return indexed_hits


def index_database_words(db_records, k):
    """
    Build a database-wide inverted k-mer index.

    This function processes every sequence in the database and records which
    database sequence IDs contain each k-mer. The resulting structure allows
    query k-mers to be looked up directly instead of scanning every database
    sequence for every query.

    Parameters
    ----------
    db_records : SequenceDatabase
        Normalized database object. Records must contain unique "id"
        values and "sequence" strings.
    k : int
        Previously validated positive k-mer size.

    Returns
    -------
    dict
        A mapping with two entries: ``kmer_index`` maps each k-mer to the set
        of database IDs that contain it, and ``sequence_lookup`` maps each ID
        to its sequence.

        Example:
        {
            "kmer_index": {
                "ATG": {"seq1", "seq2"},
                "TGC": {"seq1", "seq2"},
                "GCG": {"seq1", "seq2"},
                "CGT": {"seq1"},
                "CGA": {"seq2"}
            },
            "sequence_lookup": {
                "seq1": "ATGCGT",
                "seq2": "ATGCGA"
            }
        }

    Notes
    -----
    Each k-mer maps to a set of sequence IDs, not a list. This prevents repeated
    occurrences of the same k-mer within one database sequence from adding the
    same sequence ID multiple times to the index.

    This index is presence-based, not position-based. It records which sequences
    contain each k-mer, but not where the k-mer occurs. Positional indexing can
    be added later for seed-extension algorithms.
    """
    db_records = db_records.get_sequences()

    kmer_index = {}
    sequence_lookup = {}

    for record in db_records:
        sequence_lookup[record["id"]] = record["sequence"]

        kmers = generate_kmers(record["sequence"], k)

        for word in kmers:
            if word not in kmer_index:
                kmer_index[word] = {record["id"]}
            else:
                kmer_index[word].add(record["id"])

    return {
        "kmer_index": kmer_index,
        "sequence_lookup": sequence_lookup,
    }
