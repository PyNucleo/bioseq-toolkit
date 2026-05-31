from collections import Counter
#Dicts for transcribing template & coding strands
DNA_TEMPLATE_TO_MRNA = {"A":"U",
              "G":"C",
              "T":"A",
              "C":"G"}
DNA_CODING_TO_MRNA = {"A":"A",
              "G":"G",
              "T":"U",
              "C":"C"}

DNA_COMPLEMENT = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G"
}

def total_length(sequence):
    return len(sequence)

def base_number(sequence): 
    sequence = sequence.upper()
    return Counter(sequence)

def gc_content(sequence):

    sequence = sequence.upper()

    if not sequence:
        return 0.0
    
    counts = base_number(sequence)

    gc_fraction = (counts["G"] + counts["C"]) / total_length(sequence)
    return round(gc_fraction * 100, 2)

def reverse_complement(sequence):
    sequence = sequence.upper()
    return "".join(DNA_COMPLEMENT[base] for base in sequence)[::-1]

#Transcription:
def mrna_template(sequence):
    sequence = sequence.upper()
    return "".join(DNA_TEMPLATE_TO_MRNA[base] for base in sequence)
def mrna_coding(sequence):
    sequence = sequence.upper()
    return ''.join(DNA_CODING_TO_MRNA[base] for base in sequence)