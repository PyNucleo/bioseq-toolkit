def kmer_search(query, DB, k, threshold):

    """
     Initializes a dictionary of Database sequence and shared kmers pairs
     Generates a list of words for the query sequence of size k each, iterates over each Database sequence, gets that number of words shared between the database sequence and the query sequence;
     If the number of words shared is less than a threshold, that sequence iss eliminated, otherwise, the item of the key containing the DB sequence is updated

    After potential candidates are found, they are filtered down further based on a set ratio, where anything having words less than the max shared kmers * ratio is also eliminated.
    Result is returned as a dictionary.
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
    return set(seq[i : i + k] for i in range(len(seq) - k + 1))

def get_shared_kmers(query_words, seq, k):
    #Break down both sequences into words of size k
    
    db_seq_words = generate_kmers(seq, k)

    shared_kmers = 0

    for word in query_words:
        if word in db_seq_words:
            shared_kmers += 1

    return shared_kmers

def filter_by_relative_score(max_kmers, candidates_dict, ratio = 0.3):

    filtered = {
    seq_id: count
    for seq_id, count in candidates_dict.items()
    if count >= max_kmers * ratio
}
    return filtered


    