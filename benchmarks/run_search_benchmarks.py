from benchmarks.benchmark_search import (
    benchmark_kmer_search_run_time,
    benchmark_kmer_refinement_run_time
)
from benchmarks.benchmark_utils import (
    get_average_runtime
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
SETTINGS = [
    {"k": 3, "threshold": 1},
    {"k": 3, "threshold": 3},
    {"k": 3, "threshold": 5},
    {"k": 4, "threshold": 2},
    {"k": 4, "threshold": 3},
]

def run():
    for dataset_name, dataset_path in DATASETS.items():
        iterations = ITERATIONS[dataset_name]

        print(f"\nRuntimes for {dataset_name}:")
        print("--------------------------")

        print("\nK-mer only:")
        for setting in SETTINGS:
            k = setting["k"]
            threshold = setting["threshold"]

            runtime = get_average_runtime(
                lambda dataset, k=k, threshold=threshold: benchmark_kmer_search_run_time(
                    dataset,
                    k=k,
                    threshold=threshold,
                ),
                iterations,
                dataset_path,
            )

            print(f"k={k}, threshold={threshold}: {runtime:.6f}s")

        print("\nK-mer + SW refinement:")
        runtime = get_average_runtime(
            lambda dataset: benchmark_kmer_refinement_run_time(
                dataset,
                k=3,
                threshold=3,
                top_n_hits=10,
            ),
            iterations,
            dataset_path,
        )

        print(f"k=3, threshold=3, top_n_hits=10: {runtime:.6f}s")

if __name__ == "__main__":
    run()