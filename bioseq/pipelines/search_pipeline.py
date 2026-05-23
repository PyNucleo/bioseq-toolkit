from database.database_utils import normalize_database
from database.sequence_database import SequenceDatabase

from bioseq.search.kmer_search import kmer_search
from bioseq.search.similarity_search import rank_database_sequences


def search(query, database=None, k=3, threshold=1, top_n_hits=10, refinement=False):
    
    """
    Required Parameters: Query sequence and DataBase
    Optional Parameters: Word size (k), minimum threshold for shared words between query and DB sequence to consider the Database sequence worthy of being analyzed, maximum number of desired homologs/hits, and a boolean refinement option to be implemented further.

    Steps:

    1: Input DB is normalized from the beginning for smooth processing across the rest of the pipeline processes

    2: A k-mer search is performed on the query and the "sequence" strings in the normalized DB data structure as Biological sequences, under the assumption that related sequences must share common words (aka contiguous stretches of residues of size k)
       This returns a dictionary of "sequence" : "shared kmers" containing potential homologs (Based on the word size and minimum threshold of shared kmers).

    3: The potential homologs/hits are then stored in a list for simple processing, which is then turned into a DataBase object from the SequenceDatabase class.
       Once again, the consistency of using an Object to represent a collection/Database is maintained for consistent processing logic across all methods invoked.

    4: We pass the potential hits Database into a "rank_database_sequences" method, rearranging the sequences in a descending order based on the most optimal alignment score of each.
       This result is returned in the form of a list of tuples, with the 0th index containing "sequence" with a corresponding score at the 1th index.

    5: Finally, we trim down this list to only include up to the first top n hits only, mainting the order of each (sequence, score) tuple within the list.

    """
    db = normalize_database(database)

    potential_hits = kmer_search(query, db, k, threshold) #Dictionary of sequence : shared_kmer pairs 

    candidate_sequences = [key for key in potential_hits]
    candidate_sequences = SequenceDatabase(candidate_sequences)

    ranked_hits = rank_database_sequences(query, candidate_sequences) #Ranks based on scores; list of tuples: [(sequence, score)]



    top_hits = ranked_hits[:top_n_hits]

    return top_hits

