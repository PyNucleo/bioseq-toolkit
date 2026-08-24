# Stage A checkpoint

Recorded on 2026-08-23 before final evidence generation. The generated smoke
JSON was removed from the final artifact set after these probe results and the
smoke procedure were captured here and in deterministic tests. Smoke evidence
was never included in study summaries.

## Workload 2 — exhaustive Smith-Waterman

- Callable: `benchmarks.benchmark_alignment.benchmark_smith_waterman_scores_run_time`
- Dataset: `data/benchmark_sequences/astral_100.fasta` (100 records, 14,274 residues)
- Query: established benchmark protein query (151 residues)
- Scoring: simple match `2`, mismatch `-1`, linear gap `-2`; score-only Smith-Waterman
- One-operation probe: 1.313564 s manual; 1.234972 s pyperf retained smoke value
- pyperf smoke calibration: 1 loop, 1 worker run, 1 retained value
- Manual sample count: 50
- Estimated one manual execution: 65.7 s
- Estimated one default pyperf execution: about 120 s
- Estimated five-manual total: 328 s (5.5 min)
- Estimated five-pyperf total: 600 s (10 min)
- Estimated full workload total: about 15.5 min

## Workload 3 — k-mer filtering plus refinement

- Callable: `benchmarks.benchmark_search.benchmark_kmer_refinement_run_time`
- Dataset: `data/benchmark_sequences/astral_1000.fasta` (1,000 records, 142,363 residues)
- Query: established benchmark protein query (151 residues)
- `k`: 3
- Absolute threshold: 3
- `top_n_hits`: 10
- Refinement/scoring: enabled; simple match `1`, mismatch `-1`, linear gap `-2`
- Fixed relative candidate filter: unchanged at 0.3
- One-operation probe: 0.192885 s manual; 0.246563 s pyperf retained smoke value
- pyperf smoke calibration: 1 loop, 1 worker run, 1 retained value
- Manual sample count: 25
- Estimated one manual execution: 4.8 s
- Estimated one default pyperf execution: about 40 s
- Estimated five-manual total: 24 s
- Estimated five-pyperf total: 200 s (3.3 min)
- Estimated full workload total: about 3.7 min

## Decision

Stage B is practical. Both workloads are bound to existing callables, no
production behavior change is required, pyperf worker argument propagation and
manager-only manual execution were observed, artifact paths reject overwrite,
and the projected combined Stage-B duration is about 19–20 minutes.
