import json
import csv
import hashlib
import math
from copy import deepcopy
from pathlib import Path

import pytest

from benchmarks.pyperf_spike.analyze_spike import (
    calculate_statistics,
    evidence_rows,
    load_manual_raw,
    load_pyperf_raw,
    validate_pyperf_evidence,
)
from benchmarks.pyperf_spike.run_spike import (
    artifact_path,
    get_workload,
    refuse_existing,
    should_run_manual,
    write_json_exclusive,
)


def test_known_and_unknown_workload_selection():
    workload = get_workload("kmer_refinement")
    assert workload.parameters["top_n_hits"] == 10

    with pytest.raises(ValueError, match="unknown workload"):
        get_workload("not-a-workload")


def test_artifact_paths_separate_workload_method_and_execution():
    spike_root = Path("C:/spike-test-root")
    first = artifact_path(
        "kmer_refinement", "perf_counter", "run_01", spike_root=spike_root
    )
    second = artifact_path(
        "kmer_refinement", "pyperf", "run_01", spike_root=spike_root
    )
    third = artifact_path(
        "kmer_refinement", "perf_counter", "run_02", spike_root=spike_root
    )
    other = artifact_path(
        "exhaustive_smith_waterman",
        "perf_counter",
        "run_01",
        spike_root=spike_root,
    )

    assert len({first, second, third, other}) == 4
    assert first.name == "run_01.json"
    assert first.parent.name == "perf_counter"


def test_artifact_paths_reject_nonstudy_execution_ids():
    spike_root = Path("C:/spike-test-root")
    with pytest.raises(ValueError, match="run_01 through run_05"):
        artifact_path(
            "kmer_refinement", "perf_counter", "run_06", spike_root=spike_root
        )
    smoke = artifact_path(
        "kmer_refinement",
        "perf_counter",
        "run_00",
        smoke=True,
        spike_root=spike_root,
    )
    assert "smoke" in smoke.parts


def test_exclusive_write_and_preflight_refuse_overwrite():
    path = Path("benchmarks/pyperf_spike/smoke/test_guard/evidence.json")
    path.unlink(missing_ok=True)
    try:
        write_json_exclusive(path, {"value": 1})
        assert json.loads(path.read_text(encoding="utf-8")) == {"value": 1}

        with pytest.raises(FileExistsError, match="refusing to overwrite"):
            refuse_existing([path])
        with pytest.raises(FileExistsError):
            write_json_exclusive(path, {"value": 2})
        assert json.loads(path.read_text(encoding="utf-8")) == {"value": 1}
    finally:
        path.unlink(missing_ok=True)


def test_statistics_contract_uses_sample_standard_deviation():
    stats = calculate_statistics([1.0, 2.0, 3.0, 4.0])
    assert stats["minimum"] == 1.0
    assert stats["maximum"] == 4.0
    assert stats["median"] == 2.5
    assert stats["MAD"] == 1.0
    assert stats["mean"] == 2.5
    assert stats["sample_std"] == pytest.approx(math.sqrt(5 / 3))
    assert stats["CV_percent"] == pytest.approx(math.sqrt(5 / 3) / 2.5 * 100)


def test_manual_raw_summary_row():
    path = Path("benchmarks/pyperf_spike/smoke/test_manual_row.json")
    path.unlink(missing_ok=True)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(
                {
                    "workload": "kmer_refinement",
                    "execution_id": "run_01",
                    "durations_seconds": [1.0, 2.0, 3.0],
                }
            ),
            encoding="utf-8",
        )
        row = load_manual_raw(path)
        assert row["sample_count"] == 3
        assert row["mean"] == 2.0
        assert row["evidence_kind"] == "raw_observations"
    finally:
        path.unlink(missing_ok=True)


def test_manual_path_is_manager_only():
    assert should_run_manual("both", is_worker=False)
    assert should_run_manual("manual", is_worker=False)
    assert not should_run_manual("both", is_worker=True)
    assert not should_run_manual("pyperf", is_worker=False)


def test_historical_evidence_matches_recorded_hashes_and_sizes():
    manifest = Path(
        "benchmarks/pyperf_spike/workloads/non_refined_search/"
        "historical_evidence_sha256.csv"
    )
    with manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 11
    for row in rows:
        path = Path(row["current_path"])
        assert path.stat().st_size == int(row["bytes"])
        assert hashlib.sha256(path.read_bytes()).hexdigest().upper() == row["sha256"]


def test_study_discovery_excludes_smoke_and_preserves_historical_manual_evidence():
    for workload in (
        "non_refined_search",
        "exhaustive_smith_waterman",
        "kmer_refinement",
    ):
        rows = evidence_rows(workload)
        assert len(rows) == 10
        assert {row["execution_id"] for row in rows} == {
            "run_01",
            "run_02",
            "run_03",
            "run_04",
            "run_05",
        }
        assert all(row["evidence_kind"] != "smoke" for row in rows)

    historical = evidence_rows("non_refined_search")
    manual = [row for row in historical if row["method"] == "perf_counter"]
    assert len(manual) == 5
    assert all(row["evidence_kind"] == "historical_summary_only" for row in manual)
    assert all(row["sample_count"] == 240 for row in manual)


def test_valid_replacement_pyperf_run_has_expected_retained_structure():
    path = Path(
        "benchmarks/pyperf_spike/workloads/kmer_refinement/raw/pyperf/run_03.json"
    )
    row = load_pyperf_raw(path, "kmer_refinement")
    assert row["execution_id"] == "run_03"
    assert row["sample_count"] == 60
    assert row["evidence_kind"] == "native_pyperf_json"


def test_pyperf_integrity_rejects_chronology_contradiction_not_fast_values():
    path = Path(
        "benchmarks/pyperf_spike/workloads/kmer_refinement/raw/pyperf/run_03.json"
    )
    payload = json.loads(path.read_text(encoding="utf-8"))
    validate_pyperf_evidence(payload, path, "kmer_refinement")

    contradictory = deepcopy(payload)
    contradictory["benchmarks"][0]["runs"][2]["metadata"]["duration"] = 382.226
    with pytest.raises(ValueError, match="duration contradicts"):
        validate_pyperf_evidence(contradictory, path, "kmer_refinement")


def test_pyperf_integrity_rejects_incomplete_retained_value_structure():
    path = Path(
        "benchmarks/pyperf_spike/workloads/kmer_refinement/raw/pyperf/run_03.json"
    )
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["benchmarks"][0]["runs"][1]["values"].pop()
    with pytest.raises(ValueError, match="must contain 3 values"):
        validate_pyperf_evidence(payload, path, "kmer_refinement")
