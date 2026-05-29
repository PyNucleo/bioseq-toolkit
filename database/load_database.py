from bioseq.fasta_io import read_fasta_records
from database.sequence_database import SequenceDatabase

def create_database(db):

    records = read_fasta_records(db)

    return SequenceDatabase(records)
    