from .pipelines.translation_pipeline import process_fasta_sequences
from .search.similarity_search import rank_by_shared_kmers
from .translation import translate_sequence
from .alignment.needleman_wunsch import global_alignment


def main():
    
    action = int(input("What would you like to perform?\n"
                       "1 for similarity search\n"
                       "2 for FASTA translation pipeline\n"
                       "> "
    ))
    
    if action == 1:
        query = "ATACAGAT"
        results = rank_by_shared_kmers(query, "GATTACA")
        print(results)
    
    elif action == 2:
        file_path = input("Insert your file path: ")
        results = process_fasta_sequences(file_path)
        print(results)

if __name__ == "__main__":
    main()