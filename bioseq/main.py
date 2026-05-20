import pandas as pd

from .fasta_io import read_fasta
from .validators import validate_dna
from .sequence_utils import mrna_template, mrna_coding, total_length, gc_content, reverse_complement, base_number
from .translation import translate_dna
from .alignment.needleman_wunsch import global_alignment
from .search.similarity_search import wtvr
def main():
    
    action = int(input("What would you like to perform: \n1 for similarity search,"))
    
    if action == 1:
        wtvr("GATTACA")
        
        return
        
    FILE = input("Insert your file path: ")
    file_record = read_fasta(FILE)
    sequences = []
    mrna_seqs = []
    
    s1 = "AGT"
    s2 = "GTT"

    final = {"Sequence": None,
             "Length" : None,
             "Base counts" : None,
             "GC Content" : None,
             "Transcribed Strand" : None,
             "Amino acid chain" : None
             }
    for result in file_record:
        if result["sequence"]:
            seq = result["sequence"]
            if validate_dna(seq):
                sequences.append(seq)
                mrna_seqs.append(mrna_template(seq))
            else:
                continue
        else:
            continue
    final["Sequence"] = sequences
    final["Length"] = [total_length(s) for s in sequences]
    final["Base counts"] = [base_number(s) for s in sequences]
    final["GC Content"] = [gc_content(s) for s in sequences]
    final["Transcribed Strand"] = [mrna for mrna in mrna_seqs]
    final["Amino acid chain"] = [translate_dna(s) for s in mrna_seqs]
    #print(len(sequences))

    #df = pd.DataFrame(final)
    #df.to_csv( FILE + "test.csv", index = False)
    print(global_alignment(s1, s2))

if __name__ == "__main__":
    main()