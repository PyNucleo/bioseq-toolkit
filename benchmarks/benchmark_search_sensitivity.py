import csv
import time
from collections import Counter

from bioseq.alignment.smith_waterman import get_best_scores
from bioseq.fasta_io import read_fasta_records
from bioseq.search.kmer_search import kmer_search
from bioseq.search.similarity_search import rank_by_shared_kmers
from database.sequence_database import SequenceDatabase
from benchmarks.benchmark_utils import build_sw_reference_groups, csv_loader

dataset_10000 = "data/benchmark_sequences/astral_10000.fasta"
sw_reference_groups = None

query = (
    "ttltnpqkaairsswskfmdngvsngqgfymdlfkahpetltpfkslfggltlaqlqdnpkmkaqslvfcngmssfvdhlddndmlvvliqkmaklhnnrgirasdlrtaydilihymedhnhmvggakdawevfvgficktlgdymkels"
)
K_VALUES = [3]
THRESHOLDS = [1]
# K_VALUES = [2, 3, 4, 5]
# THRESHOLDS = [1, 2, 3, 5, 10]

MATRIX = "BLOSUM62"
GAP_PENALTY = -4

USE_CACHED_SW_REFERENCE = True
SW_REFERENCE_CSV = "benchmarks/cached_exhaustive_sw_results/astral_10000_BLOSUM62_gap-4.csv"
CACHED_SW_10000_RUNTIME = 639.6809130999973

OUTPUT_CSV = "benchmarks/search_sensitivity_results_astral_10000_k3_threshold1.csv"

def get_hit_ids(results):
    return [hit["id"] for hit in results]


def sort_results(results, score_key):
    return sorted(results, key=lambda hit: hit[score_key], reverse=True)


def get_tier_name(seq_id):
    for tier_name, ids in sw_reference_groups.items():
        if seq_id in ids:
            return tier_name

    return "unclassified"


def count_tiers(results):
    tier_counts = Counter()

    for hit in results:
        tier_counts[get_tier_name(hit["id"])] += 1

    return tier_counts


def exhaustive_sw_search(query, records, matrix=None, match=2, mismatch=-1, gap_penalty=-4):
    start = time.perf_counter()

    results = []

    for record in records:
        score, best_positions = get_best_scores(
            query,
            record["sequence"],
            gap_penalty=gap_penalty,
            matrix=matrix,
            match=match,
            mismatch=mismatch,
        )

        results.append(
            {
                "id": record["id"],
                "sequence": record["sequence"],
                "score": score,
                "best_positions": best_positions,
            }
        )

    runtime = time.perf_counter() - start

    return sort_results(results, "score"), runtime


def kmer_only_search(query, records, k=3, threshold=1):
    start = time.perf_counter()

    db = SequenceDatabase(records)

    candidates = kmer_search(
        query=query,
        db=db,
        k=k,
        threshold=threshold,
    )

    ranked_hits = rank_by_shared_kmers(query, candidates)

    runtime = time.perf_counter() - start

    return ranked_hits, runtime


def refined_kmer_search(
    query,
    kmer_hits,
    matrix=None,
    match=2,
    mismatch=-1,
    gap_penalty=-4,
):
    start = time.perf_counter()

    refined_hits = []

    for hit in kmer_hits:
        score, best_positions = get_best_scores(
            query,
            hit["sequence"],
            gap_penalty=gap_penalty,
            matrix=matrix,
            match=match,
            mismatch=mismatch,
        )

        refined_hit = hit.copy()
        refined_hit["sw_score"] = score
        refined_hit["best_positions"] = best_positions

        refined_hits.append(refined_hit)

    refined_hits = sort_results(refined_hits, "sw_score")

    runtime = time.perf_counter() - start

    return refined_hits, runtime


def recall_count(reference_results, method_results, n):
    reference_ids = set(get_hit_ids(reference_results[:n]))
    method_ids = set(get_hit_ids(method_results[:n]))

    return len(reference_ids & method_ids)


def high_confidence_recall_count(method_results):
    method_ids = set(get_hit_ids(method_results))

    tier_1_ids = set(sw_reference_groups["tier_1_exact_or_near_exact"])
    tier_2_ids = set(sw_reference_groups["tier_2_high_scoring"])
    high_confidence_ids = tier_1_ids | tier_2_ids

    recovered = method_ids & high_confidence_ids

    return len(recovered), len(high_confidence_ids)


def background_returned_count(method_results):
    method_ids = set(get_hit_ids(method_results))
    background_ids = set(sw_reference_groups["tier_4_background"])

    return len(method_ids & background_ids)


def make_sweep_row(
    k,
    threshold,
    database_size,
    exhaustive_sw_results,
    sw_runtime,
):
    kmer_results, kmer_runtime = kmer_only_search(
        query=query,
        records=records,
        k=k,
        threshold=threshold,
    )

    refined_results, refined_runtime = refined_kmer_search(
        query=query,
        kmer_hits=kmer_results,
        matrix=MATRIX,
        gap_penalty=GAP_PENALTY,
    )

    refined_total_runtime = kmer_runtime + refined_runtime

    kmer_high_recovered, high_total = high_confidence_recall_count(kmer_results)
    refined_high_recovered, _ = high_confidence_recall_count(refined_results)

    candidate_count = len(kmer_results)
    candidate_fraction = candidate_count / database_size if database_size > 0 else 0

    kmer_speedup = sw_runtime / kmer_runtime if kmer_runtime > 0 else float("inf")
    refined_speedup = (
        sw_runtime / refined_total_runtime
        if refined_total_runtime > 0
        else float("inf")
    )

    kmer_tiers = count_tiers(kmer_results)
    refined_tiers = count_tiers(refined_results)

    return {
        "k": k,
        "threshold": threshold,
        "candidate_count": candidate_count,
        "candidate_fraction": candidate_fraction,
        "kmer_recall_5": f"{recall_count(exhaustive_sw_results, kmer_results, 5)}/5",
        "kmer_recall_10": f"{recall_count(exhaustive_sw_results, kmer_results, 10)}/10",
        "kmer_recall_20": f"{recall_count(exhaustive_sw_results, kmer_results, 20)}/20",
        "kmer_high_confidence_recall": f"{kmer_high_recovered}/{high_total}",
        "refined_recall_5": f"{recall_count(exhaustive_sw_results, refined_results, 5)}/5",
        "refined_recall_10": f"{recall_count(exhaustive_sw_results, refined_results, 10)}/10",
        "refined_recall_20": f"{recall_count(exhaustive_sw_results, refined_results, 20)}/20",
        "refined_high_confidence_recall": f"{refined_high_recovered}/{high_total}",
        "background_returned": background_returned_count(kmer_results),
        "kmer_tier_1": kmer_tiers["tier_1_exact_or_near_exact"],
        "kmer_tier_2": kmer_tiers["tier_2_high_scoring"],
        "kmer_tier_3": kmer_tiers["tier_3_low_moderate"],
        "kmer_tier_4": kmer_tiers["tier_4_background"],
        "refined_tier_1": refined_tiers["tier_1_exact_or_near_exact"],
        "refined_tier_2": refined_tiers["tier_2_high_scoring"],
        "refined_tier_3": refined_tiers["tier_3_low_moderate"],
        "refined_tier_4": refined_tiers["tier_4_background"],
        "sw_runtime": sw_runtime,
        "kmer_runtime": kmer_runtime,
        "refined_total_runtime": refined_total_runtime,
        "kmer_speedup": kmer_speedup,
        "refined_speedup": refined_speedup,
    }


def run_benchmark():
    global records

    records = read_fasta_records(dataset_10000)
    database_size = len(records)

    if USE_CACHED_SW_REFERENCE:
        print("Loading cached exhaustive SW results...")
        exhaustive_sw_results = csv_loader(SW_REFERENCE_CSV)
        sw_runtime = CACHED_SW_10000_RUNTIME
    else:
        print("Running exhaustive Smith-Waterman once...")
        
        exhaustive_sw_results, sw_runtime = exhaustive_sw_search(
            query=query,
            records=records,
            matrix=MATRIX,
            gap_penalty=GAP_PENALTY,
        )
    
    global sw_reference_groups

    sw_reference_groups = build_sw_reference_groups(exhaustive_sw_results)

    fieldnames = [
        "k",
        "threshold",
        "candidate_count",
        "candidate_fraction",
        "kmer_recall_5",
        "kmer_recall_10",
        "kmer_recall_20",
        "kmer_high_confidence_recall",
        "refined_recall_5",
        "refined_recall_10",
        "refined_recall_20",
        "refined_high_confidence_recall",
        "background_returned",
        "kmer_tier_1",
        "kmer_tier_2",
        "kmer_tier_3",
        "kmer_tier_4",
        "refined_tier_1",
        "refined_tier_2",
        "refined_tier_3",
        "refined_tier_4",
        "sw_runtime",
        "kmer_runtime",
        "refined_total_runtime",
        "kmer_speedup",
        "refined_speedup",
    ]

    rows = []

    for k in K_VALUES:
        for threshold in THRESHOLDS:
            row = make_sweep_row(
                k=k,
                threshold=threshold,
                database_size=database_size,
                exhaustive_sw_results=exhaustive_sw_results,
                sw_runtime=sw_runtime,
            )

            rows.append(row)

            print(
                f"done: k={k}, threshold={threshold}, "
                f"candidates={row['candidate_count']}, "
                f"kmer_recall@20={row['kmer_recall_20']}, "
                f"refined_recall@20={row['refined_recall_20']}"
            )

    with open(OUTPUT_CSV, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print()
    print(f"Saved sweep results to: {OUTPUT_CSV}")


if __name__ == "__main__":
    run_benchmark()