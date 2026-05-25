from bioseq.fasta_io import read_fasta_sequences_only
from database.sequence_database import SequenceDatabase

def create_database(db):

    sequences = read_fasta_sequences_only(db)

    db = SequenceDatabase()

    return SequenceDatabase(sequences)
    