# pyperf evaluation spike

## Purpose and decision

This completed research spike evaluates whether pyperf provides enough
methodological value to become Bioseq Toolkit's standard timing layer for
deliberate benchmark evidence. It compares the project's complete simple
`time.perf_counter()` procedure with pyperf's complete benchmark procedure; it
is not merely a comparison between timer functions.

The project-level decision is:

- use `time.perf_counter()` for quick development timing, reconnaissance, and
  rough implementation before/after checks;
- use pyperf for deliberate reproducible benchmark evidence, formal runtime
  characterization, and the future regular-versus-indexed benchmark.

The choice depends on benchmark purpose, not algorithm type. pyperf supplies
warmups, automatic calibration, retained measurements, worker processes,
runtime/environment metadata, diagnostics, a standardized measurement
structure, and raw native evidence that can be reanalyzed later.

## Installation and reproduction

From the repository root, install the project with its development dependencies:

```powershell
python -m pip install -e ".[dev]"
```

The `dev` extra reproducibly declares `pyperf==2.10.0`, the version used by the
preserved study evidence. Confirm the active environment when needed:

```powershell
python -c "import pyperf; print(pyperf.__version__)"
```

Reproduce the derived study results from the committed raw evidence without
rerunning benchmark measurements:

```powershell
python -m benchmarks.pyperf_spike.analyze_spike
python -m pytest tests/benchmarks/test_pyperf_spike.py -q
```

To execute a bounded smoke validation of the harness on a clean checkout, use
the reserved `run_00` path, which is excluded from final study summaries:

```powershell
python benchmarks/pyperf_spike/run_spike.py --workload kmer_refinement --execution run_00 --method both --smoke --debug-single-value
```

The smoke command writes generated JSON under `benchmarks/pyperf_spike/smoke/`
and refuses to overwrite an existing `run_00`. The completed `run_01` through
`run_05` study evidence should not be rerun or deleted merely to reproduce the
analysis; the harness likewise refuses to overwrite those preserved artifacts.

## Workloads

The study covers exactly three runtime regimes.

| Workload | Regime | Timed operation | Dataset | Scale | Query | Manual observations per execution |
| --- | --- | --- | --- | --- | ---: | ---: |
| `non_refined_search` | fast/noise-sensitive | public `search(query, database_path)` with refinement disabled | `data/benchmark_sequences/astral_1000.fasta` | 1,000 records; 142,363 residues | 142 residues | 240 |
| `exhaustive_smith_waterman` | expensive | `benchmarks.benchmark_alignment.benchmark_smith_waterman_scores_run_time` | `data/benchmark_sequences/astral_100.fasta` | 100 records; 14,274 residues | 151 residues | 50 |
| `kmer_refinement` | intermediate/representative | `benchmarks.benchmark_search.benchmark_kmer_refinement_run_time`, which wraps public `search()` with refinement enabled | `data/benchmark_sequences/astral_1000.fasta` | 1,000 records; 142,363 residues | 151 residues | 25 |

The exhaustive workload includes FASTA/database preparation and exhaustive
Smith-Waterman **score computation** across all 100 records. It does not perform
full traceback/alignment reconstruction. Its scoring is match `2`, mismatch
`-1`, `matrix=None`, and linear gap penalty `-2`.

The refinement workload uses `k=3`, threshold `3`, `top_n_hits=10`, refinement
enabled, and the unchanged fixed relative candidate filter `0.3`. Its scoring
is match `1`, mismatch `-1`, `matrix=None`, and linear gap penalty `-2`. FASTA
loading remains inside both current timed benchmark callables.

## Independent-execution methodology

Each workload is represented by five independent executions per timing
methodology (`run_01` through `run_05`). The manual sample budgets above reflect
operation cost and established benchmark iteration conventions; they are not
pyperf loop counts.

For workloads 2 and 3, every pyperf execution used:

- 20 retained worker runs;
- one warmup per retained worker;
- three retained values per worker, or 60 values per execution;
- automatic calibration;
- a naturally selected loop count of one.

The historical workload-1 pyperf executions have the same worker/warmup/value
shape and naturally calibrated to four loops. Its five native pyperf JSON files
retain all 60 normalized values per execution.

The harness parses pyperf's real process arguments, propagates custom workload
arguments to workers, and guards manual sampling to the manager process. Raw
paths use exclusive creation, so existing evidence cannot be overwritten.

## Historical workload-1 limitation

Workload 1 was completed before the later automated spike work and was not
rerun. Its five native pyperf files survive. The original 240 individual manual
observations per execution do not survive; only five historical manual
execution summaries remain. They are stored as
`historical_manual_summaries/` and explicitly analyzed as
`historical_summary_only`. Missing observations were neither invented nor
regenerated.

The historical harness, analysis script, original millisecond summary, and
SHA-256 inventory remain alongside that workload as provenance. The redundant
root-level legacy `summary.csv` was removed; the classified historical copy is
`workloads/non_refined_search/historical_summary_ms.csv`.

## Invalid and replacement workload-3 run_03

One original `kmer_refinement` pyperf `run_03` was excluded because the evidence
contained an independently identifiable metadata-integrity contradiction: a
worker-reported duration of hundreds of seconds was incompatible with the
sequential worker `date` and `uptime` progression. It was rejected for
structural evidence inconsistency, not because its value was statistically
inconvenient. Its raw JSON is absent and its measurements and wall-clock
duration contribute to no final result.

The current
`workloads/kmer_refinement/raw/pyperf/run_03.json` is the valid replacement. It
used the unchanged workload, benchmark parameters, harness, and normal pyperf
calibration, which again selected one loop. It contains 20 retained workers and
60 retained values. It is retained even though its mean is faster than several
other independent launches; no statistical-outlier deletion rule is applied.

The analyzer now performs a bounded structural integrity check before loading
native pyperf evidence. It requires the study's explicit calibration/worker/
warmup/value shape, finite positive retained timings, required chronology
metadata, and a worker duration consistent with sequential `date` and `uptime`
progression. It raises an error and never deletes evidence. It does not reject
high CVs, slow or fast runs, or means far from other launches.

## Final results

The following values are regenerated from the current valid raw evidence.
Execution means are the mean of the five independent execution means.
Across-execution CV is calculated over those five means.

| Workload | Manual execution mean | pyperf execution mean | Manual avg within-execution CV | pyperf avg within-execution CV | Relative within-CV reduction | Manual across-execution CV | pyperf across-execution CV |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Non-refined search | 40.664 ms | 38.767 ms | 12.008% | 3.576% | 70.2% | 1.246% | 0.888% |
| Exhaustive Smith-Waterman score computation | 1.2709 s | 1.2458 s | 3.715% | 3.351% | 9.8% | 2.088% | 5.361% |
| K-mer filtering + Smith-Waterman refinement | 189.292 ms | 162.948 ms | 6.638% | 4.161% | 37.3% | 0.490% | 14.957% |

For the fast workload, pyperf produced substantially tighter within-execution
distributions and somewhat stronger repeatability across independent launches.
For exhaustive Smith-Waterman score computation, the within-execution
improvement was modest and across-launch repeatability was worse. For the
combined refinement workload, the average within-execution distribution was
tighter while independent pyperf launches varied substantially more than the
manual execution means.

Differences between manual and pyperf mean timings do **not** show that pyperf
made an algorithm faster. The methodologies have different measurement
structures. This spike evaluates within-execution distributions,
across-launch repeatability, warmups, calibration, isolation, metadata,
diagnostics, and raw-evidence auditability—not an algorithm-speed causal effect.

## Statistics and reproducibility

Per-execution summaries use seconds and report minimum, maximum, median,
median absolute deviation (`median(abs(x_i - median(x)))`), mean, sample
standard deviation (the `n - 1` denominator), coefficient of variation
(`sample_std / mean * 100`), and `maximum / median` as a descriptive slow-tail
indicator. Relative within-execution CV reduction is
`((manual_CV - pyperf_CV) / manual_CV) * 100`.

Regenerate derived summaries without running benchmark workloads:

```powershell
python -m benchmarks.pyperf_spike.analyze_spike
```

The analyzer reads only the canonical workload evidence directories. It does
not discover `smoke/`, root exploratory JSON, or derived/legacy summaries as
study evidence. Generated smoke JSON was removed after the probe methodology
and results were preserved in `STAGE_A_CHECKPOINT.md` and tests.

`execution_durations.csv` records outer-command wall-clock measurements when
they were captured. The rejected original workload-3 `run_03` duration was
removed. The exact outer-command duration of the valid replacement was not
recorded, so that cell is blank with an explicit note rather than a fabricated
value. Internal pyperf worker duration metadata is not the same measurement and
was not substituted.

## Artifact layout

```text
benchmarks/pyperf_spike/
  README.md
  STAGE_A_CHECKPOINT.md
  run_spike.py
  analyze_spike.py
  execution_durations.csv
  overall_summary.csv
  workloads/
    non_refined_search/
      historical_harness.py
      historical_analysis.py
      historical_evidence_sha256.csv
      historical_summary_ms.csv
      historical_manual_summaries/run_NN.json
      raw/pyperf/run_NN.json
      summary.csv
    exhaustive_smith_waterman/
      raw/perf_counter/run_NN.json
      raw/pyperf/run_NN.json
      summary.csv
    kmer_refinement/
      raw/perf_counter/run_NN.json
      raw/pyperf/run_NN.json
      summary.csv
```

`STAGE_A_CHECKPOINT.md` is retained because it records workload-selection
rationale, probe durations, sample-budget estimates, and the practicality
decision that preceded final evidence generation.

## Limitations and future implication

pyperf improves and characterizes measurement behavior within an execution;
it does not eliminate environmental variation between independent benchmark
launches. The replacement `run_03` is direct evidence that a structurally valid
launch can differ materially from other sessions.

Future regular-versus-indexed timing should therefore use closely paired,
controlled comparisons rather than comparing unrelated sessions. That future
design should retain a correctness/equivalence gate, deliberate order handling,
database size and query count as dimensions, separate index-build and reuse
costs, raw machine-readable evidence, dataset/query hashes, and versioned
configuration. None of that future benchmark is implemented by this spike.

This study is methodological only. It does not establish biological usefulness,
homology-detection quality, statistical biological significance, BLAST
equivalence, production readiness, or indexed-search superiority.
