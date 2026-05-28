from bioseq.fasta_io import read_fasta_sequences_only, read_fasta_records
from database.sequence_database import SequenceDatabase

FILE = "C:/Users/Admin/Documents/TemporaryHoldings(CanBeDeleted)/astral-scopedom-seqres-gd-all-2.08-stable.fa.txt"

def dataset_to_chunks(dataset, *chunk_lengths):
    
    dataset_sequences = read_fasta_records(FILE)

    for chunk_size in chunk_lengths:
        temp_chunk = dataset_sequences[:chunk_size]

        with open("data/benchmark_sequences/astral_" + str(chunk_size)  + ".txt", "w") as file:
            
            for record in temp_chunk:
                file.write(record["header"] + "\n")
                file.write(record["sequence"] + "\n")


dataset_to_chunks(FILE, 10, 100, 1000, 10000)

