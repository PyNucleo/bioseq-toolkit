from .pipelines.translation_pipeline import process_fasta_sequences
from .pipelines.search_pipeline import search
from .search.similarity_search import rank_by_shared_kmers
from .translation import translate_sequence
from .alignment.needleman_wunsch import global_alignment


def main():
    
    print("Run a small demo search example.")
    
    results = search(
        query="ATGCG",
        database=[
            "ATGCGA",
            "ATGCGA",
            "GGGGGG",
        ],
        k=3,
        threshold=1,
        refinement=True,
    )

    for hit in results:
        print(hit)

if __name__ == "__main__":
    main()