"""Bounded perf_counter-versus-pyperf study harness.

The harness deliberately supports only the three pyperf-spike workloads.  The
first workload is historical and cannot be rerun through this script.
"""

from __future__ import annotations

import json
import os
import platform
import re
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

import pyperf


SPIKE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = SPIKE_ROOT.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmarks.benchmark_alignment import (  # noqa: E402
    benchmark_smith_waterman_scores_run_time,
)
from benchmarks.benchmark_search import (  # noqa: E402
    benchmark_kmer_refinement_run_time,
)


HISTORICAL_QUERY = (
    "sltsadkshvrsiwskaggsaeeigaealgrmlesfpntktyfdhyadlsvssaqvhthgkkiidalttav"
    "nhidditgalsslstlhaqtlrvdpanfkilshtilvvlalyfpadftpevhlacdkflanvshaladnyr"
)
BENCHMARK_QUERY = (
    "ttltnpqkaairsswskfmdngvsngqgfymdlfkahpetltpfkslfggltlaqlqdnpkmkaqslvfc"
    "ngmssfvdhlddndmlvvliqkmaklhnnrgirasdlrtaydilihymedhnhmvggakdawevfvgfick"
    "tlgdymkels"
)
EXECUTION_RE = re.compile(r"run_(\d{2})\Z")


@dataclass(frozen=True)
class Workload:
    name: str
    label: str
    runtime_regime: str
    callable_name: str
    dataset: str
    dataset_records: int
    dataset_residues: int
    query: str
    manual_samples: int
    scoring: str
    parameters: dict[str, object]
    historical: bool = False

    @property
    def query_length(self) -> int:
        return len(self.query)

    @property
    def dataset_path(self) -> Path:
        return REPO_ROOT / self.dataset


WORKLOADS = {
    "non_refined_search": Workload(
        name="non_refined_search",
        label="Non-refined public search",
        runtime_regime="fast",
        callable_name="bioseq.pipelines.search_pipeline.search",
        dataset="data/benchmark_sequences/astral_1000.fasta",
        dataset_records=1000,
        dataset_residues=142363,
        query=HISTORICAL_QUERY,
        manual_samples=240,
        scoring="inactive because refinement=False",
        parameters={"refinement": False},
        historical=True,
    ),
    "exhaustive_smith_waterman": Workload(
        name="exhaustive_smith_waterman",
        label="Exhaustive Smith-Waterman score search",
        runtime_regime="expensive",
        callable_name=(
            "benchmarks.benchmark_alignment."
            "benchmark_smith_waterman_scores_run_time"
        ),
        dataset="data/benchmark_sequences/astral_100.fasta",
        dataset_records=100,
        dataset_residues=14274,
        query=BENCHMARK_QUERY,
        manual_samples=50,
        scoring="simple match=2, mismatch=-1, linear gap=-2",
        parameters={
            "match": 2,
            "mismatch": -1,
            "gap_penalty": -2,
            "matrix": None,
            "alignment": "Smith-Waterman score only",
        },
    ),
    "kmer_refinement": Workload(
        name="kmer_refinement",
        label="K-mer filtering plus Smith-Waterman refinement",
        runtime_regime="intermediate",
        callable_name=(
            "benchmarks.benchmark_search.benchmark_kmer_refinement_run_time"
        ),
        dataset="data/benchmark_sequences/astral_1000.fasta",
        dataset_records=1000,
        dataset_residues=142363,
        query=BENCHMARK_QUERY,
        manual_samples=25,
        scoring="simple match=1, mismatch=-1, linear gap=-2",
        parameters={
            "k": 3,
            "threshold": 3,
            "top_n_hits": 10,
            "refinement": True,
            "relative_candidate_filter": 0.3,
            "match_score": 1,
            "mismatch_score": -1,
            "gap_penalty": -2,
            "matrix": None,
        },
    ),
}


def get_workload(name: str) -> Workload:
    """Resolve a known workload or fail with the available names."""
    try:
        return WORKLOADS[name]
    except KeyError as exc:
        choices = ", ".join(sorted(WORKLOADS))
        raise ValueError(f"unknown workload {name!r}; choose one of: {choices}") from exc


def validate_execution_id(execution_id: str, *, smoke: bool) -> None:
    match = EXECUTION_RE.fullmatch(execution_id)
    if match is None:
        raise ValueError("execution ID must use run_NN format")
    number = int(match.group(1))
    if smoke:
        if number != 0:
            raise ValueError("smoke executions must use run_00")
    elif number not in range(1, 6):
        raise ValueError("study executions must be run_01 through run_05")


def artifact_path(
    workload: str,
    method: str,
    execution_id: str,
    *,
    smoke: bool = False,
    spike_root: Path = SPIKE_ROOT,
) -> Path:
    """Return a unique semantic artifact path without creating it."""
    get_workload(workload)
    if method not in {"perf_counter", "pyperf"}:
        raise ValueError("method must be 'perf_counter' or 'pyperf'")
    validate_execution_id(execution_id, smoke=smoke)
    if smoke:
        return spike_root / "smoke" / workload / method / f"{execution_id}.json"
    return (
        spike_root
        / "workloads"
        / workload
        / "raw"
        / method
        / f"{execution_id}.json"
    )


def refuse_existing(paths: list[Path]) -> None:
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError(
            "refusing to overwrite existing evidence: " + ", ".join(existing)
        )


def write_json_exclusive(path: Path, payload: dict[str, object]) -> None:
    """Write JSON using exclusive creation so races cannot overwrite evidence."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def should_run_manual(method: str, *, is_worker: bool) -> bool:
    return not is_worker and method in {"manual", "both"}


def operation_for(workload: Workload) -> Callable[[], object]:
    """Bind the existing project callable without reimplementing its work."""
    dataset = str(workload.dataset_path)
    if workload.name == "exhaustive_smith_waterman":
        return lambda: benchmark_smith_waterman_scores_run_time(dataset)
    if workload.name == "kmer_refinement":
        return lambda: benchmark_kmer_refinement_run_time(
            dataset,
            k=3,
            threshold=3,
            top_n_hits=10,
            refinement=True,
        )
    raise ValueError("the historical non_refined_search workload is not rerunnable")


def git_revision() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def workload_metadata(workload: Workload) -> dict[str, object]:
    data = asdict(workload)
    data["query_length"] = workload.query_length
    data["dataset_path"] = str(workload.dataset_path)
    return data


def run_manual(
    workload: Workload,
    operation: Callable[[], object],
    execution_id: str,
    output: Path,
    *,
    smoke: bool,
) -> None:
    sample_count = 1 if smoke else workload.manual_samples
    started = datetime.now(timezone.utc)
    wall_start = time.perf_counter()
    durations = []
    for _ in range(sample_count):
        start = time.perf_counter()
        operation()
        durations.append(time.perf_counter() - start)
    wall_duration = time.perf_counter() - wall_start
    completed = datetime.now(timezone.utc)

    payload = {
        "schema_version": 1,
        "workload": workload.name,
        "execution_id": execution_id,
        "timing_method": "perf_counter",
        "evidence_kind": "smoke" if smoke else "study",
        "started_utc": started.isoformat(),
        "completed_utc": completed.isoformat(),
        "wall_clock_duration_seconds": wall_duration,
        "environment": {
            "git_revision": git_revision(),
            "python_version": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "platform": platform.platform(),
            "pyperf_version": pyperf.__version__,
            "process_id": os.getpid(),
        },
        "benchmark": workload_metadata(workload),
        "sample_count": sample_count,
        "unit": "second",
        "durations_seconds": durations,
    }
    write_json_exclusive(output, payload)
    print(f"manual evidence: {output} ({sample_count} observations)")


def _add_worker_arguments(cmd: list[str], args: object) -> None:
    cmd.extend(
        [
            "--workload",
            args.workload,
            "--execution",
            args.execution,
            "--method",
            args.method,
        ]
    )
    if args.smoke:
        cmd.append("--smoke")


def main() -> int:
    runner = pyperf.Runner(add_cmdline_args=_add_worker_arguments)
    runner.argparser.add_argument("--workload", required=True, choices=WORKLOADS)
    runner.argparser.add_argument("--execution", required=True)
    runner.argparser.add_argument(
        "--method", choices=("manual", "pyperf", "both"), default="both"
    )
    runner.argparser.add_argument(
        "--smoke",
        action="store_true",
        help="write run_00 outside final evidence and use one manual observation",
    )

    # Parse the real process arguments.  pyperf workers therefore retain their
    # --worker, --worker-task, --loops, warmup, and value controls.
    args = runner.parse_args()
    workload = get_workload(args.workload)
    validate_execution_id(args.execution, smoke=args.smoke)
    if workload.historical:
        raise ValueError("non_refined_search is preserved historical evidence")
    if args.output or args.append:
        raise ValueError("output paths are managed by the spike harness")

    operation = operation_for(workload)
    manual_path = artifact_path(
        workload.name, "perf_counter", args.execution, smoke=args.smoke
    )
    pyperf_path = artifact_path(
        workload.name, "pyperf", args.execution, smoke=args.smoke
    )

    if not args.worker:
        requested = []
        if args.method in {"manual", "both"}:
            requested.append(manual_path)
        if args.method in {"pyperf", "both"}:
            requested.append(pyperf_path)
        refuse_existing(requested)

    if should_run_manual(args.method, is_worker=args.worker):
        run_manual(
            workload,
            operation,
            args.execution,
            manual_path,
            smoke=args.smoke,
        )

    if args.method == "manual":
        if args.worker:
            raise RuntimeError("a pyperf worker cannot run the manual-only method")
        return 0

    runner.metadata.update(
        {
            "workload": workload.name,
            "execution_id": args.execution,
            "timing_methodology": "pyperf",
            "evidence_kind": "smoke" if args.smoke else "study",
            "git_revision": git_revision(),
            "dataset": workload.dataset,
            "dataset_records": workload.dataset_records,
            "dataset_residues": workload.dataset_residues,
            "query_length": workload.query_length,
            "callable": workload.callable_name,
            "scoring": workload.scoring,
            "parameters_json": json.dumps(workload.parameters, sort_keys=True),
        }
    )
    if not args.worker:
        pyperf_path.parent.mkdir(parents=True, exist_ok=True)
        args.output = str(pyperf_path)

    runner.bench_func(f"{workload.name}_pyperf", operation)
    if not args.worker:
        print(f"pyperf evidence: {pyperf_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
