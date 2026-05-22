def kmer_search(queue, DB, k, threshold):

    max_shared_kmers = 0
    sequences_and_shared_kmers = {x : 0 for x in DB}

    for db_seq in DB:
        temp_shared_kmers = get_shared_kmers(queue, db_seq, k)

        if temp_shared_kmers < threshold:
            sequences_and_shared_kmers.pop(db_seq)
        
        else:
            if temp_shared_kmers > max_shared_kmers:
                max_shared_kmers = temp_shared_kmers

            sequences_and_shared_kmers[db_seq] = temp_shared_kmers

    return filter_candidates_ratio(max_shared_kmers, sequences_and_shared_kmers)


def generate_kmers(seq, k):
    return set(seq[i : i + k] for i in range(len(seq) - k + 1))

def get_shared_kmers(queue, seq, k):
    #Break down both sequences into words of size k
    queue_words = generate_kmers(queue, k)
    db_seq_words = generate_kmers(seq, k)

    shared_kmers = 0

    for word in queue_words:
        if word in db_seq_words:
            shared_kmers += 1

    return shared_kmers

def filter_candidates_ratio(max_kmers, candidates_dict, ratio = 0.3):

    filtered = {
    seq_id: count
    for seq_id, count in candidates_dict.items()
    if count >= max_kmers * ratio
}
    return filtered


    