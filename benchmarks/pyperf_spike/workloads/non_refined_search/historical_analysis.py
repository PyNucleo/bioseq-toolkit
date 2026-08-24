import csv
import json
import pandas as pd
from experiment_pyperf import get_stats
rows = []

for execution in range(1, 6):

    perf_counter_json = f"benchmarks/pyperf_spike/executions_{execution}_perf_counter.json"
    pyperf_json = f"benchmarks/pyperf_spike/executions_{execution}_pyperf.json"

    with open(perf_counter_json, 'r', encoding='utf-8') as perf_counter_file:
        perf_counter_data = json.load(perf_counter_file)
        perf_counter_data["execution"] = execution
        perf_counter_data["method"] = "perf_counter"

    with open(pyperf_json, 'r', encoding='utf-8') as pyperf_file:
        pyperf_data = json.load(pyperf_file)

        pyperf_times = []

        for run in pyperf_data["benchmarks"][0]["runs"]:
            if "values" in run:
                pyperf_times.extend(value * 1000 for value in run["values"])

        pyperf_stats = get_stats(pyperf_times)
        pyperf_stats["execution"] = execution
        pyperf_stats["method"] = "pyperf"

    rows.append(perf_counter_data)
    rows.append(pyperf_stats)

df = pd.DataFrame(rows)[
    [
        "execution",
        "method",
        "minimum",
        "maximum",
        "median",
        "MAD",
        "mean",
        "std",
        "CV",
    ]
]

df.to_csv("benchmarks/pyperf_spike/summary.csv", index=False)
