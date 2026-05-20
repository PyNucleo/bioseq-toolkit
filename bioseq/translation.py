from Bio.Seq import Seq

def translate_dna(seq, stop_at_stop = True):
    seq = Seq(seq)
    aa = seq.translate(to_stop = stop_at_stop)
    return aa
def translate_rna(seq, stop_at_stop = True):
    seq = Seq(seq)
    aa = seq.translate(to_stop = stop_at_stop)
    return aa