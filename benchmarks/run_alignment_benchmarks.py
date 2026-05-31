from benchmarks.benchmark_alignment import (
    benchmark_smith_waterman_scores_run_time,
    benchmark_smith_waterman_alignment_run_time,
)
from benchmarks.benchmark_utils import (
    get_average_runtime,
    get_total_dataset_residues,
)


DATASETS = {
    "astral_10": "data/benchmark_sequences/astral_10.fasta",
    "astral_100": "data/benchmark_sequences/astral_100.fasta",
    "astral_1000": "data/benchmark_sequences/astral_1000.fasta",
    "astral_10000": "data/benchmark_sequences/astral_10000.fasta",
}

ITERATIONS = {
    "astral_10": 100,
    "astral_100": 50,
    "astral_1000": 25,
    "astral_10000": 10,
}


def run():
    print("dataset,mode,iterations,total_residues,average_runtime_seconds")

    for dataset_name, dataset_path in DATASETS.items():
        iterations = ITERATIONS[dataset_name]
        residues = get_total_dataset_residues(dataset_path)

        score_runtime = get_average_runtime(
            benchmark_smith_waterman_scores_run_time,
            iterations,
            dataset_path,
        )

        alignment_runtime = get_average_runtime(
            benchmark_smith_waterman_alignment_run_time,
            iterations,
            dataset_path,
        )

        print(f"{dataset_name},score_only,{iterations},{residues},{score_runtime:.6f}")
        print(f"{dataset_name},alignment_reconstruction,{iterations},{residues},{alignment_runtime:.6f}")


if __name__ == "__main__":
    run()