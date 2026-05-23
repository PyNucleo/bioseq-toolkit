from Bio.Seq import Seq

def translate_sequence(seq, stop_at_stop = True):
    return Seq(seq).translate(to_stop = stop_at_stop) 
