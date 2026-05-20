from ..fasta_io import read_fasta
from ..validators import valid_dna
from ..sequence_utils import(
   mrna_template,
    total_length,
    gc_content,
    base_number,
)

from ..translation import translate_dna


def process_fasta_sequences(file_path):
    file_records = read_fasta(file_path)


    sequences = []
    mrna_seqs = []

    for record in file_records:
        seq = record.get("sequence")

        if not seq:
            continue

        if not valid_dna(seq):
            continue

        sequences.append(seq)
        mrna_seqs.append(mrna_template(seq))
    
    final = {
        "Sequence": sequences,
        "Length": [total_length(seq) for seq in sequences],
        "Base Counts": [gc_content(seq) for seq in sequences],
        "Transcribed Strand": mrna_seqs,
        "Amino acid chain": [translate_dna(mrna) for mrna in mrna_seqs],
    }

    return final