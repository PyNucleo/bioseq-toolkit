import pyperf
import time
import json
import statistics

from bioseq.pipelines.search_pipeline import search


query = "sltsadkshvrsiwskaggsaeeigaealgrmlesfpntktyfdhyadlsvssaqvhthgkkiidalttavnhidditgalsslstlhaqtlrvdpanfkilshtilvvlalyfpadftpevhlacdkflanvshaladnyr"
database = "C:/Users/Admin/bioseq-toolkit/data/benchmark_sequences/astral_1000.fasta"


def get_stats(recorded_times):
    minimum = min(recorded_times)
    maximum = max(recorded_times)
    median = statistics.median(recorded_times)
    mad = statistics.median(abs(x - median) for x in recorded_times)
    mean = statistics.mean(recorded_times)
    std = statistics.stdev(recorded_times)
    cv = (std / mean) * 100

    return {
        "minimum": minimum,
        "maximum": maximum,
        "median": median,
        "MAD": mad,
        "mean": mean,
        "std": std,
        "CV": cv,
    }
def get_perf_counter_stats(loop_iterations):
    recorded_times = []

    for _ in range(loop_iterations):
        start = time.perf_counter()
        search(query, database)
        end = time.perf_counter()

        recorded_times.append((end - start) * 1000)

    return get_stats(recorded_times)

if __name__ == "__main__":
    execution = 5

    runner = pyperf.Runner()

    # Read the real arguments supplied to this process.
    # This also allows pyperf workers to receive their --worker,
    # --loops, etc. arguments correctly.
    runner.parse_args()

    # Only the main pyperf process should run our manual benchmark.
    if not runner.args.worker:
        perf_counter_json_path = (
            f"C:/Users/Admin/bioseq-toolkit/"
            f"benchmarks/pyperf_spike/"
            f"executions_{execution}_perf_counter.json"
        )

        perf_counter_metadata = get_perf_counter_stats(240)

        with open(perf_counter_json_path, "w", encoding="utf-8") as f:
            json.dump(perf_counter_metadata, f, indent=4)

    # Both the main process and pyperf workers must reach this.
    runner.bench_func(
        "perf counter vs pyperf benchmark",
        search,
        query,
        database,
    )