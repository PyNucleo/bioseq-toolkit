from database.database_utils import normalize_database
from database.sequence_database import SequenceDatabase

from bioseq.search.kmer_search import kmer_search
from bioseq.search.similarity_search import rank_database_sequences


def search(query, database=None, k=3, threshold=1, top_n_hits=10, refinement=False):
    
    db = normalize_database(database)

    potential_hits = kmer_search(query, db, k, threshold) #Dictionary of sequence : shared_kmer pairs 

    potential_sequences = [key for key in potential_hits]
    potential_sequences = SequenceDatabase(potential_sequences)

    ranked_hits = rank_database_sequences(query, potential_sequences) #Ranks based on scores; sequence : score pairs



    top_hits = ranked_hits[:top_n_hits]

    return top_hits

