import pandas as pd

from .pipelines.translation_pipeline import process_fasta_sequences
from .search.similarity_search import rank_database_sequences
from .translation import translate_dna
from .alignment.needleman_wunsch import global_alignment
from .search.similarity_search import rank_database_sequences
def main():
    
    action = int(input("What would you like to perform?\n"
                       "1 for similarity search\n"
                       "2 for FASTA translation pipeline\n"
                       "> "
    ))
    
    if action == 1:
        results = rank_database_sequences("GATTACA")
        print(results)
    
    elif action == 2:
        file_path = input("Insert your file path: ")
        results = process_fasta_sequences(file_path)
        print(results)

if __name__ == "__main__":
    main()