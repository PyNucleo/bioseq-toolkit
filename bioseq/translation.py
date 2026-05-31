from Bio.Seq import Seq

def translate_sequence(seq, stop_at_stop=True, trim_partial=True):
    if trim_partial:
        usable_length = len(seq) - (len(seq) % 3)
        seq = seq[:usable_length]

    return Seq(seq).translate(to_stop=stop_at_stop)
