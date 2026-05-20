from collections import Counter
#Dicts for transcribing template & coding strands
DICT_CONVERTTemplate = {"A":"U",
              "G":"C",
              "T":"A",
              "C":"G"}
DICT_CONVERTCoding = {"A":"A",
              "G":"G",
              "T":"U",
              "C":"C"}

def total_length(DNA):
    return len(DNA)

def base_number(DNA): 
    return Counter(DNA)

def gc_content(DNA):
    S = base_number(DNA)
    C_Content = S["C"]
    G_Content = S["G"]
    GC_Hundo = (C_Content + G_Content) / total_length(DNA)
    return f"{GC_Hundo * 100:.2f}"

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
        else:
            E+="G"
    return E

#Transcription:
def mrna_template(DNA_Sequence):
    return ''.join(DICT_CONVERTTemplate[base] for base in DNA_Sequence)
def mrna_coding(DNA_Sequence):
    return ''.join(DICT_CONVERTCoding[base] for base in DNA_Sequence)