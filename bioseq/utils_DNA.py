from collections import Counter
from Bio import SeqIO

import pandas as pd



Valid_DNA = ["A", "T", "G", "C"]
Valid_mRNA = ["A", "U", "G", "C"]
DICT_CONVERTTemplate = {"A":"U",
              "G":"C",
              "T":"A",
              "C":"G"}
DICT_CONVERTCoding = {"A":"A",
              "G":"A",
              "T":"U",
              "C":"C"}
CODON_DICT = {}
Phe_Key=["UUU","UUC"]
Leu_Key=["UUA","UUG","CUU","CUC","CUA","CUG"]
Ile_Key=["AUU","AUC","AUA"]
Val_Key=["GUU","GUC","GUA","GUG"]
Ser_Key=["UCU","UCC","UCA","UCG", "AGU","AGC"]
Pro_Key=["CCU","CCC","CCA","CCG"]
Thr_Key=["ACU","ACC","ACA","ACG"]
Ala_Key=["GCU","GCC","GCA","GCG"]
Tyr_Key=["UAU","UAC"]
His_Key=["CAU","CAC"]
Gln_Key=["CAA","CAG"]
Asn_Key=["AAU","AAC"]
Lys_Key=["AAA","AAG"]
Asp_Key=["GAU","GAC"]
Glu_Key=["GAA","GAG"]
Cys_Key=["UGU","UGC"]
Trp_Key=["UGG"]
Arg_Key=["CGU","CGC","CGA","CGG", "AGA","AGG"]
Gly_Key=["GGU","GGC","GGA","GGG"]
STOP_Key=["UAA","UAG","UGA"]
Met_Key=["AUG"]
CODON_DICT.update({Phe:"Phe" for Phe in Phe_Key})
CODON_DICT.update({Ile:"Ile" for Ile in Ile_Key})
CODON_DICT.update({Val:"Val" for Val in Val_Key})
CODON_DICT.update({Ser:"Ser" for Ser in Ser_Key})         
CODON_DICT.update({Pro:"Pro" for Pro in Pro_Key})            
CODON_DICT.update({Leu:"Leu" for Leu in Leu_Key})
CODON_DICT.update({Thr:"Thr" for Thr in Thr_Key})         
CODON_DICT.update({Ala:"Ala" for Ala in Ala_Key})
CODON_DICT.update({Tyr:"Tyr" for Tyr in Tyr_Key})
CODON_DICT.update({His:"His" for His in His_Key})
CODON_DICT.update({Gln:"Gln" for Gln in Gln_Key})
CODON_DICT.update({Asn:"Asn" for Asn in Asn_Key})
CODON_DICT.update({Lys:"Lys" for Lys in Lys_Key})            
CODON_DICT.update({Asp:"Asp" for Asp in Asp_Key})
CODON_DICT.update({Glu:"Glu" for Glu in Glu_Key})
CODON_DICT.update({Cys:"Cys" for Cys in Cys_Key})
CODON_DICT.update({Trp:"Trp" for Trp in Trp_Key})
CODON_DICT.update({Arg:"Arg" for Arg in Arg_Key})
CODON_DICT.update({Gly:"Gly" for Gly in Gly_Key})
CODON_DICT.update({STOP:"STOP" for STOP in STOP_Key})
CODON_DICT.update({Met:"Met" for Met in Met_Key})

def invalid_DNA(DNA_Sequence):
        return any(i not in Valid_DNA for i in set(DNA_Sequence))

def invalid_mRNA(mRNA_Sequence):
        return any(i not in Valid_mRNA for i in set(mRNA_Sequence))
    
def total_length(DNA):
    return len(DNA)

def base_number(DNA):
    Base_Counts = Counter(DNA) 
    return Base_Counts

def gc_content_per(DNA):
    S = base_number(DNA)
    C_Content=S["C"]
    G_Content=S["G"]
    GC_Hundo = ((C_Content+G_Content)/total_length(DNA))
    return GC_Hundo * 100

def reverse_complement(DNA):
    DNA = DNA[::-1]
    E = ''
    for comp in range(0, len(DNA)):
        if DNA[comp]=="A":
            E+="T"
        elif DNA[comp]=="T":
            E+="A"
        elif DNA[comp]=="G":
            E+="C"
        elif DNA[comp]=="C":
            E+="G"
    return E


def mrna_template(DNA_Sequence):
    return ''.join(DICT_CONVERTTemplate[base] for base in DNA_Sequence)
def mrna_coding(DNA_Sequence):
    return ''.join(DICT_CONVERTCoding[base] for base in DNA_Sequence)

def protein(DNA_Sequence):
    AA_Seq = []
    l = len(DNA_Sequence)
    mRNA = list(mrna_coding(DNA_Sequence))
    for i in range(0, l - (l % 3), 3):
        c_out = mRNA[i:i + 3]
        s=''.join(c_out)
        AA_Seq.append(CODON_DICT[s])

        if CODON_DICT[s] == "STOP":
            return AA_Seq
    return AA_Seq

def analyze(DNA_Sequence, n):
    all_results = []
    results = {}
    if n == "1":       
        results["Sequence"] = DNA_Sequence
        results['Length'] = total_length(DNA_Sequence)
        results['Base counts'] = base_number(DNA_Sequence)
        results['GC content'] = gc_content_per(DNA_Sequence)
        #results['Reverse complement'] = reverse_complement(DNA_Sequence)
        results['Transcription'] = mrna_coding(DNA_Sequence)
        results['Amino acid chain'] = protein(DNA_Sequence)
    elif n == "2":
        results["Sequence"] = DNA_Sequence
        results['Length'] = total_length(DNA_Sequence)
        results['Base counts'] = base_number(DNA_Sequence)
        results['GC content'] = gc_content_per(DNA_Sequence)
        results['Reverse complement'] = reverse_complement(DNA_Sequence)
        results['Transcription'] = mrna_template(DNA_Sequence)
        results['Amino acid chain'] = protein(DNA_Sequence)
        
    else:
        results["Sequence"] = DNA_Sequence
        results['Length'] = total_length(DNA_Sequence)
        results['Base counts'] = base_number(DNA_Sequence)
        results['GC content'] = gc_content_per(DNA_Sequence)
        results['Amino acid chain'] = protein(DNA_Sequence)
    return results
def main():
    sequence_type = input("What is your sequence type? Is it: \n" + 
                          "1: DNA Coding strand \n" + 
                          "2: DNA Template strand \n" +
                          "3: mRNA \n")
    FILE = input("File name: ").strip()
    final = []
    sequence = ''
    sequences = [] 
    with open(FILE, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    sequences.append(sequence)
                    sequence = ''
            else:
                sequence += line
        sequences.append(sequence)
    for DNA_Sequence in sequences:
        if not invalid_DNA(sequence):
            final.append(analyze(DNA_Sequence, sequence_type))
    db = pd.DataFrame(final)
    db.to_csv("test1.csv", index = False)
main()