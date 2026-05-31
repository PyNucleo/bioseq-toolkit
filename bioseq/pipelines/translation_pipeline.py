from ..fasta_io import read_fasta_records
from ..validators import valid_dna
from ..translation import translate_sequence
from ..sequence_utils import(
   mrna_template,
    total_length,
    gc_content,
    base_number,
)


def process_fasta_sequences(file_path):
    file_records = read_fasta_records(file_path)

    sequences = []
    mrna_seqs = []

    for record in file_records:
        seq = record["sequence"].upper()

        if not seq:
            continue

        if not valid_dna(seq):
            continue

        sequences.append(seq)
        mrna_seqs.append(mrna_template(seq))
    
    final = {
        "sequence": sequences,
        "length": [total_length(seq) for seq in sequences],
        "base_counts": [gc_content(seq) for seq in sequences],
        "transcribed_strand": mrna_seqs,
        "amino_acid_chain": [translate_sequence(mrna) for mrna in mrna_seqs],
    }

    return final