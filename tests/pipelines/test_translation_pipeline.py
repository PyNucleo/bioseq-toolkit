from bioseq.pipelines.translation_pipeline import process_fasta_sequences
from bioseq.validators import valid_dna

def test_lower_and_upper_case_dna(tmp_path):

    fasta = tmp_path / "tests.fasta"

    fasta.write_text(
        ">seq1\n"
        "atgc\n"
        ">seq2\n"
        "ac\n"
        ">seq3\n"
        "GC\n"
        ">seq4\n"
        "c\n"
        ">seq5\n"
        "actcgcagtagca\n"
        ">seq6\n"
        "ATGGCCATT\n"
    )

    result =  {
        "sequence": ["ATGC","AC","GC","C","ACTCGCAGTAGCA", "ATGGCCATT"],
        "length": [4, 2, 2, 1, 13, 9],
        "gc_content": [50, 50, 100, 100, 53.85, 44.44],
        "transcribed_strand": ["UACG", "UG","CG","G","UGAGCGUCAUCGU", "UACCGGUAA"],
        "amino_acid_chain": ["Y", "", "", "", "", "YR"],
    }

    assert process_fasta_sequences(fasta) == result