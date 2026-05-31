from bioseq.alignment.smith_waterman import local_alignment, get_best_scores
from database.database_utils import normalize_database
from bioseq.fasta_io import read_fasta_sequences_only, read_fasta_records

import time

def benchmark_smith_waterman_alignment_run_search(FILE):
    
    sw_alignments = []

    query = "ttltnpqkaairsswskfmdngvsngqgfymdlfkahpetltpfkslfggltlaqlqdnpkmkaqslvfcngmssfvdhlddndmlvvliqkmaklhnnrgirasdlrtaydilihymedhnhmvggakdawevfvgficktlgdymkels"
    db = normalize_database(read_fasta_sequences_only(FILE))

    db_sequences = db.get_sequences()

    for seq2 in db_sequences:
            sw_alignments.append({
            "sequence":seq2,
            "alignments":local_alignment(
                            query,
                            seq2["sequence"],
                            match=2,
                            mismatch=-1,
                            gap_penalty=-2,
                            return_all=False,
                            structured=True,
                        )
        })

    return sw_alignments

def benchmark_smith_waterman_alignment_run_time(FILE):
    
    start = time.perf_counter()


    query = "ttltnpqkaairsswskfmdngvsngqgfymdlfkahpetltpfkslfggltlaqlqdnpkmkaqslvfcngmssfvdhlddndmlvvliqkmaklhnnrgirasdlrtaydilihymedhnhmvggakdawevfvgficktlgdymkels"
    db = normalize_database(read_fasta_sequences_only(FILE))

    db_sequences = db.get_sequences()

    for seq2 in db_sequences:
        local_alignment(
            query,
            seq2["sequence"],
            match=2,
            mismatch=-1,
            gap_penalty=-2,
            return_all=False,
            structured=True,
        )
        

    end = time.perf_counter()

    interval = end - start

    return interval

def benchmark_smith_waterman_scores_run_results(FILE):
    start = time.perf_counter()

    all_scores = []

    query = "ttltnpqkaairsswskfmdngvsngqgfymdlfkahpetltpfkslfggltlaqlqdnpkmkaqslvfcngmssfvdhlddndmlvvliqkmaklhnnrgirasdlrtaydilihymedhnhmvggakdawevfvgficktlgdymkels"
    db = normalize_database(read_fasta_sequences_only(FILE))

    db = db.get_sequences()

    for seq in db:
        score, _ = get_best_scores(
                       query,
                       seq["sequence"],
                       gap_penalty=-2,
                       match=2,
                       mismatch=-1,
                    )
        all_scores.append({
             "sequence": seq["sequence"],
             "sw_score": score
        })

    end = time.perf_counter()

    interval = end - start

    return interval

def benchmark_smith_waterman_scores_run_time(FILE):
    start = time.perf_counter()

    query = "ttltnpqkaairsswskfmdngvsngqgfymdlfkahpetltpfkslfggltlaqlqdnpkmkaqslvfcngmssfvdhlddndmlvvliqkmaklhnnrgirasdlrtaydilihymedhnhmvggakdawevfvgficktlgdymkels"
    
    db = normalize_database(read_fasta_sequences_only(FILE))
    db = db.get_sequences()

    for seq in db:
        get_best_scores(
            query,
            seq["sequence"],
            gap_penalty=-2,
            match=2,
            mismatch=-1,
        )

    end = time.perf_counter()

    interval = end - start

    return interval