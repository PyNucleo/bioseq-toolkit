from bioseq.fasta_io import read_fasta_sequences_only
from database.sequence_database import SequenceDatabase

def create_database(DB_CHOICE):

    sequences = read_fasta_sequences_only(DB_CHOICE)

    db = SequenceDatabase()

    db.add_sequence(seq for seq in sequences)

    return db
    