import csv
from bioseq.fasta_io import read_fasta_records

def get_average_runtime(benchmark_method, iterations, dataset):
    """
    Run a benchmark function multiple times and return
    the average runtime.

    Parameters
    ----------
    benchmark_method : callable
        Benchmark function to execute.
    iterations : int
        Number of repeated runs.
    dataset : str
        Dataset path passed into benchmark_method.

    Returns
    -------
    float
        Mean runtime across all benchmark iterations.
    """

    total_runtime = 0

    for _ in range(iterations):
        runtime = benchmark_method(dataset)
        total_runtime += runtime

    return total_runtime / iterations

def get_total_dataset_residues(dataset):
    dataset_records = read_fasta_records(dataset)

    temp_s = 0

    for seq in dataset_records:
        temp_s += len(seq["sequence"])
    
    return temp_s

def build_sw_reference_groups(sw_results, top_tier_score=None, second_tier_min_ratio=0.65):
    """
    Build SW reference groups automatically from exhaustive SW results.

    Parameters
    ----------
    sw_results : list[dict]
        Exhaustive SW results sorted or unsorted.
        Each dict must contain:
        - "id"
        - "score"

    top_tier_score : float | None
        Score used for Tier 1. If None, the maximum SW score is used.

    second_tier_min_ratio : float
        Minimum score ratio relative to max_score for Tier 2.
        Example: 0.65 means Tier 2 includes hits with score >= 65% of max_score,
        excluding Tier 1.

    Returns
    -------
    dict
        Dictionary of tier names mapped to sequence IDs.
    """
    if not sw_results:
        return {
            "tier_1_exact_or_near_exact": [],
            "tier_2_high_scoring": [],
            "tier_3_low_moderate": [],
            "tier_4_background": [],
        }

    sorted_results = sorted(sw_results, key=lambda hit: hit["score"], reverse=True)

    max_score = sorted_results[0]["score"]

    if top_tier_score is None:
        top_tier_score = max_score

    tier_1 = []
    tier_2 = []
    tier_3 = []
    tier_4 = []

    tier_2_cutoff = max_score * second_tier_min_ratio

    for hit in sorted_results:
        seq_id = hit["id"]
        score = hit["score"]

        if score == top_tier_score:
            tier_1.append(seq_id)
        elif score >= tier_2_cutoff:
            tier_2.append(seq_id)
        elif score >= 100:
            tier_3.append(seq_id)
        else:
            tier_4.append(seq_id)

    return {
        "tier_1_exact_or_near_exact": tier_1,
        "tier_2_high_scoring": tier_2,
        "tier_3_low_moderate": tier_3,
        "tier_4_background": tier_4,
    }

def csv_loader(csv_path):

    results = []

    with open(csv_path, newline="") as file:
        reader = csv.DictReader(file)

        for row in reader:
            row["score"] = float(row["score"])
            row["length"] = int(row["length"])

            results.append(row)

    return results