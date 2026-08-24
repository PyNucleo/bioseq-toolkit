"""Derive reproducible spike summaries from preserved evidence."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from datetime import datetime
from pathlib import Path

import pyperf

from benchmarks.pyperf_spike.run_spike import SPIKE_ROOT, WORKLOADS


SUMMARY_FIELDS = [
    "workload",
    "method",
    "execution_id",
    "evidence_kind",
    "sample_count",
    "unit",
    "minimum",
    "maximum",
    "median",
    "MAD",
    "mean",
    "sample_std",
    "CV_percent",
    "max_over_median",
]

EXPECTED_PYPERF_VALUE_RUNS = 20
EXPECTED_PYPERF_WARMUPS_PER_RUN = 1
EXPECTED_PYPERF_VALUES_PER_RUN = 3
CHRONOLOGY_TOLERANCE_SECONDS = 0.01


def calculate_statistics(values: list[float]) -> dict[str, float | int | str]:
    if not values:
        raise ValueError("at least one retained timing value is required")
    median = statistics.median(values)
    mean = statistics.mean(values)
    sample_std = statistics.stdev(values) if len(values) > 1 else 0.0
    return {
        "sample_count": len(values),
        "unit": "second",
        "minimum": min(values),
        "maximum": max(values),
        "median": median,
        "MAD": statistics.median(abs(value - median) for value in values),
        "mean": mean,
        "sample_std": sample_std,
        "CV_percent": (sample_std / mean) * 100 if mean else 0.0,
        "max_over_median": max(values) / median if median else math.inf,
    }


def load_manual_raw(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    stats = calculate_statistics(payload["durations_seconds"])
    return {
        "workload": payload["workload"],
        "method": "perf_counter",
        "execution_id": payload["execution_id"],
        "evidence_kind": "raw_observations",
        **stats,
    }


def load_historical_manual_summary(path: Path) -> dict[str, object]:
    """Normalize the historical millisecond summary without inventing samples."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    divisor = 1000.0
    median = payload["median"] / divisor
    maximum = payload["maximum"] / divisor
    return {
        "workload": "non_refined_search",
        "method": "perf_counter",
        "execution_id": path.stem,
        "evidence_kind": "historical_summary_only",
        "sample_count": 240,
        "unit": "second",
        "minimum": payload["minimum"] / divisor,
        "maximum": maximum,
        "median": median,
        "MAD": payload["MAD"] / divisor,
        "mean": payload["mean"] / divisor,
        "sample_std": payload["std"] / divisor,
        "CV_percent": payload["CV"],
        "max_over_median": maximum / median,
    }


def validate_pyperf_evidence(
    payload: dict[str, object], path: Path, workload: str
) -> None:
    """Reject structurally inconsistent study evidence, never statistical outliers."""
    benchmarks = payload.get("benchmarks")
    if not isinstance(benchmarks, list) or len(benchmarks) != 1:
        raise ValueError(f"{path}: expected exactly one pyperf benchmark")

    benchmark = benchmarks[0]
    if not isinstance(benchmark, dict):
        raise ValueError(f"{path}: malformed pyperf benchmark payload")
    runs = benchmark.get("runs")
    if not isinstance(runs, list):
        raise ValueError(f"{path}: missing pyperf runs")

    suite_metadata = payload.get("metadata")
    if not isinstance(suite_metadata, dict):
        raise ValueError(f"{path}: missing pyperf suite metadata")
    loops = suite_metadata.get("loops")
    if isinstance(loops, bool) or not isinstance(loops, int) or loops <= 0:
        raise ValueError(f"{path}: pyperf loops must be a positive integer")

    expected_metadata = {
        "workload": workload,
        "execution_id": path.stem,
        "evidence_kind": "study",
    }
    for key, expected in expected_metadata.items():
        if key in suite_metadata and suite_metadata[key] != expected:
            raise ValueError(
                f"{path}: pyperf metadata {key!r} is "
                f"{suite_metadata[key]!r}, expected {expected!r}"
            )

    value_run_count = 0
    previous_date: datetime | None = None
    previous_uptime: float | None = None
    for index, run in enumerate(runs, start=1):
        if not isinstance(run, dict):
            raise ValueError(f"{path}: pyperf run {index} is malformed")
        metadata = run.get("metadata")
        if not isinstance(metadata, dict):
            raise ValueError(f"{path}: pyperf run {index} is missing metadata")

        try:
            current_date = datetime.fromisoformat(str(metadata["date"]))
            duration = float(metadata["duration"])
            current_uptime = float(metadata["uptime"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"{path}: pyperf run {index} has incomplete chronology metadata"
            ) from exc
        if not all(
            math.isfinite(value) and value >= 0
            for value in (duration, current_uptime)
        ):
            raise ValueError(
                f"{path}: pyperf run {index} has invalid chronology metadata"
            )

        if previous_date is not None and previous_uptime is not None:
            date_delta = (current_date - previous_date).total_seconds()
            uptime_delta = current_uptime - previous_uptime
            if date_delta < 0 or uptime_delta < 0:
                raise ValueError(
                    f"{path}: pyperf run {index} chronology moves backwards"
                )
            if (
                duration > date_delta + CHRONOLOGY_TOLERANCE_SECONDS
                and duration > uptime_delta + CHRONOLOGY_TOLERANCE_SECONDS
            ):
                raise ValueError(
                    f"{path}: pyperf run {index} duration contradicts both "
                    "sequential date and uptime progression"
                )
        previous_date = current_date
        previous_uptime = current_uptime

        if "values" not in run:
            continue
        value_run_count += 1
        values = run["values"]
        warmups = run.get("warmups")
        if not isinstance(values, list) or len(values) != EXPECTED_PYPERF_VALUES_PER_RUN:
            raise ValueError(
                f"{path}: pyperf retained run {value_run_count} must contain "
                f"{EXPECTED_PYPERF_VALUES_PER_RUN} values"
            )
        if not isinstance(warmups, list) or len(warmups) != EXPECTED_PYPERF_WARMUPS_PER_RUN:
            raise ValueError(
                f"{path}: pyperf retained run {value_run_count} must contain "
                f"{EXPECTED_PYPERF_WARMUPS_PER_RUN} warmup"
            )
        try:
            numeric_values = [float(value) for value in values]
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"{path}: pyperf retained run {value_run_count} has non-numeric values"
            ) from exc
        if not all(math.isfinite(value) and value > 0 for value in numeric_values):
            raise ValueError(
                f"{path}: pyperf retained run {value_run_count} values must be "
                "finite and positive"
            )

    if value_run_count != EXPECTED_PYPERF_VALUE_RUNS:
        raise ValueError(
            f"{path}: expected {EXPECTED_PYPERF_VALUE_RUNS} retained pyperf runs, "
            f"found {value_run_count}"
        )
    calibration_runs = len(runs) - value_run_count
    if calibration_runs != 1:
        raise ValueError(
            f"{path}: expected exactly one pyperf calibration run, "
            f"found {calibration_runs}"
        )


def load_pyperf_raw(path: Path, workload: str) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    validate_pyperf_evidence(payload, path, workload)
    benchmark = pyperf.Benchmark.load(str(path))
    values = list(benchmark.get_values())
    return {
        "workload": workload,
        "method": "pyperf",
        "execution_id": path.stem,
        "evidence_kind": "native_pyperf_json",
        **calculate_statistics(values),
    }


def evidence_rows(workload: str) -> list[dict[str, object]]:
    workload_root = SPIKE_ROOT / "workloads" / workload
    rows = []
    if workload == "non_refined_search":
        manual_dir = workload_root / "historical_manual_summaries"
        rows.extend(
            load_historical_manual_summary(path)
            for path in sorted(manual_dir.glob("run_*.json"))
        )
    else:
        manual_dir = workload_root / "raw" / "perf_counter"
        rows.extend(
            load_manual_raw(path) for path in sorted(manual_dir.glob("run_*.json"))
        )

    pyperf_dir = workload_root / "raw" / "pyperf"
    rows.extend(
        load_pyperf_raw(path, workload)
        for path in sorted(pyperf_dir.glob("run_*.json"))
    )
    return sorted(rows, key=lambda row: (str(row["execution_id"]), str(row["method"])))


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize_workload(workload: str) -> list[dict[str, object]]:
    rows = evidence_rows(workload)
    write_csv(SPIKE_ROOT / "workloads" / workload / "summary.csv", SUMMARY_FIELDS, rows)
    return rows


def _across_execution(rows: list[dict[str, object]]) -> dict[str, float]:
    means = [float(row["mean"]) for row in rows]
    sample_std = statistics.stdev(means) if len(means) > 1 else 0.0
    overall_mean = statistics.mean(means)
    return {
        "mean": overall_mean,
        "sample_std": sample_std,
        "CV_percent": (sample_std / overall_mean) * 100 if overall_mean else 0.0,
        "minimum": min(means),
        "maximum": max(means),
    }


OVERALL_FIELDS = [
    "workload",
    "runtime_regime",
    "approximate_operation_seconds",
    "manual_execution_mean_seconds",
    "pyperf_execution_mean_seconds",
    "manual_average_within_execution_CV_percent",
    "pyperf_average_within_execution_CV_percent",
    "relative_CV_reduction_percent",
    "manual_across_execution_CV_percent",
    "pyperf_across_execution_CV_percent",
    "manual_across_execution_sample_std_seconds",
    "pyperf_across_execution_sample_std_seconds",
    "manual_across_execution_minimum_mean_seconds",
    "manual_across_execution_maximum_mean_seconds",
    "pyperf_across_execution_minimum_mean_seconds",
    "pyperf_across_execution_maximum_mean_seconds",
    "manual_average_max_over_median",
    "pyperf_average_max_over_median",
    "observation",
]


def overall_row(workload: str, rows: list[dict[str, object]]) -> dict[str, object]:
    manual = [row for row in rows if row["method"] == "perf_counter"]
    perf = [row for row in rows if row["method"] == "pyperf"]
    if not manual or not perf:
        raise ValueError(f"{workload} needs both methodologies for overall analysis")
    manual_cv = statistics.mean(float(row["CV_percent"]) for row in manual)
    perf_cv = statistics.mean(float(row["CV_percent"]) for row in perf)
    relative = ((manual_cv - perf_cv) / manual_cv) * 100 if manual_cv else ""
    manual_across = _across_execution(manual)
    perf_across = _across_execution(perf)
    within = "lower" if perf_cv < manual_cv else "not lower"
    across = (
        "lower"
        if perf_across["CV_percent"] < manual_across["CV_percent"]
        else "not lower"
    )
    observation = (
        f"pyperf within-execution CV was {within}; "
        f"pyperf across-execution CV was {across}."
    )
    all_means = [float(row["mean"]) for row in rows]
    return {
        "workload": workload,
        "runtime_regime": WORKLOADS[workload].runtime_regime,
        "approximate_operation_seconds": statistics.mean(all_means),
        "manual_execution_mean_seconds": manual_across["mean"],
        "pyperf_execution_mean_seconds": perf_across["mean"],
        "manual_average_within_execution_CV_percent": manual_cv,
        "pyperf_average_within_execution_CV_percent": perf_cv,
        "relative_CV_reduction_percent": relative,
        "manual_across_execution_CV_percent": manual_across["CV_percent"],
        "pyperf_across_execution_CV_percent": perf_across["CV_percent"],
        "manual_across_execution_sample_std_seconds": manual_across["sample_std"],
        "pyperf_across_execution_sample_std_seconds": perf_across["sample_std"],
        "manual_across_execution_minimum_mean_seconds": manual_across["minimum"],
        "manual_across_execution_maximum_mean_seconds": manual_across["maximum"],
        "pyperf_across_execution_minimum_mean_seconds": perf_across["minimum"],
        "pyperf_across_execution_maximum_mean_seconds": perf_across["maximum"],
        "manual_average_max_over_median": statistics.mean(
            float(row["max_over_median"]) for row in manual
        ),
        "pyperf_average_max_over_median": statistics.mean(
            float(row["max_over_median"]) for row in perf
        ),
        "observation": observation,
    }


def generate(workloads: list[str]) -> None:
    all_rows = {workload: summarize_workload(workload) for workload in workloads}
    if set(workloads) == set(WORKLOADS):
        rows = [overall_row(name, all_rows[name]) for name in WORKLOADS]
        write_csv(SPIKE_ROOT / "overall_summary.csv", OVERALL_FIELDS, rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workload", choices=(*WORKLOADS, "all"), default="all")
    args = parser.parse_args()
    names = list(WORKLOADS) if args.workload == "all" else [args.workload]
    generate(names)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
