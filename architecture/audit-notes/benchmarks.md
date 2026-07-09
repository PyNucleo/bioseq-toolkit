# Benchmarks

## Purpose

Document current benchmark scripts, datasets, reports, and their relationship to public search/alignment behavior. Benchmark outputs are measurement artifacts for current code, inputs, hardware, and parameters; they are not automatic proof of general biological performance.

## Location

- `benchmarks/benchmark_alignment.py`
- `benchmarks/benchmark_search.py`
- `benchmarks/benchmark_search_sensitivity.py`
- `benchmarks/benchmark_utils.py`
- `benchmarks/run_alignment_benchmarks.py`
- `benchmarks/run_search_benchmarks.py`
- `benchmarks/main.py`
- `benchmarks/BENCHMARKS.md`
- `benchmarks/SEARCH_SENSITIVITY.md`
- `data/benchmark_sequences/*.fasta`
- `benchmarks/cached_exhaustive_sw_results/*.csv`
- `benchmarks/figures/*.png`

## Public entry points

- `benchmark_smith_waterman_alignment_run_time(FILE)`
- `benchmark_smith_waterman_scores_run_time(FILE)`
- `benchmark_kmer_search_run_time(file, k=3, threshold=3)`
- `benchmark_kmer_refinement_run_time(file, k=3, threshold=1, top_n_hits=10, refinement=True)`
- `run_alignment_benchmarks.run()`
- `run_search_benchmarks.run()`
- `benchmark_search_sensitivity.run_benchmark()`
- Utility helpers such as `get_average_runtime`, `get_total_dataset_residues`, `build_sw_reference_groups`, and `csv_loader`.

## Inputs

- Benchmark FASTA paths under `data/benchmark_sequences`.
- Hard-coded protein query string in benchmark scripts.
- Benchmark settings such as `k`, `threshold`, `top_n_hits`, matrix, and gap penalty.
- Optional cached exhaustive Smith-Waterman CSV for sensitivity benchmark.

## Input normalization and validation

- Alignment benchmarks read FASTA records or sequence-only lists through FASTA helpers and normalize sequence-only lists into `SequenceDatabase`.
- Search benchmarks pass file paths through `normalize_database`, so they are parsed as FASTA database inputs.
- Sensitivity benchmark reads FASTA records directly and wraps them in `SequenceDatabase` for k-mer search.
- `csv_loader` casts cached CSV `score` to `float` and `length` to `int`.
- No benchmark script validates biological representativeness of the data beyond using the current files.

## Main algorithm or workflow

Alignment runtime benchmarks:

1. Read benchmark FASTA.
2. Use one hard-coded protein query.
3. Compare full local alignment reconstruction with score-only `get_best_scores`.
4. Repeat via `get_average_runtime` in runner script.
5. Print CSV-style rows.

Search runtime benchmarks:

1. Read/normalize benchmark FASTA.
2. Measure scan-based `kmer_search` across settings.
3. Measure pipeline search with Smith-Waterman refinement on top hits.
4. Print runtime summaries.

Sensitivity benchmark:

1. Load records from `astral_10000.fasta`.
2. Load cached exhaustive Smith-Waterman results when enabled, or compute exhaustive results.
3. Build score-derived reference tiers.
4. Run k-mer-only search and refined k-mer search.
5. Compute recall counts, candidate fractions, tier counts, background returned, and speedups.
6. Write a CSV report.

## Data structures created or consumed

- Benchmark dataset records: FASTA-derived record dictionaries with `id`, metadata, and `sequence`.
- Runtime return values: `float` seconds from `time.perf_counter`.
- Exhaustive SW result dictionaries:
  - `id`: record ID.
  - `sequence`: sequence string.
  - `score`: best Smith-Waterman score.
  - `best_positions`: best matrix positions.
- Sensitivity row dictionary:
  - settings: `k`, `threshold`.
  - candidate metrics: `candidate_count`, `candidate_fraction`.
  - recall strings: `kmer_recall_5`, `kmer_recall_10`, `kmer_recall_20`, and refined equivalents.
  - tier counts: `kmer_tier_*`, `refined_tier_*`.
  - runtimes: `sw_runtime`, `kmer_runtime`, `refined_total_runtime`.
  - speedups: `kmer_speedup`, `refined_speedup`.
- Reference tier dictionary from `build_sw_reference_groups`:
  - `tier_1_exact_or_near_exact`: `list[str]` IDs.
  - `tier_2_high_scoring`: `list[str]` IDs.
  - `tier_3_low_moderate`: `list[str]` IDs.
  - `tier_4_background`: `list[str]` IDs.

## Return value

- Runtime benchmark functions return `float`.
- Some result-producing helpers return lists of dictionaries and runtime tuples.
- Runner scripts print results.
- Sensitivity benchmark writes CSV output and prints progress.

## Side effects

- Benchmark runners print to stdout.
- Sensitivity benchmark writes `benchmarks/search_sensitivity_results_astral_10000_k3_threshold1.csv`.
- Dataset chunking utility outside benchmarks writes FASTA subsets under `data/benchmark_sequences`.
- Benchmark scripts read local dataset and cached CSV files.

## Dependencies called

- `bioseq.alignment.smith_waterman.local_alignment`
- `bioseq.alignment.smith_waterman.get_best_scores`
- `bioseq.search.kmer_search.kmer_search`
- `bioseq.pipelines.search_pipeline.search`
- `bioseq.fasta_io.read_fasta_records`
- `bioseq.fasta_io.read_fasta_sequences_only`
- `database.database_utils.normalize_database`
- Python `time`, `csv`, and `collections.Counter`.

## Assumptions and limitations

- Observed: benchmark query is hard-coded in multiple files.
- Observed: current reports state hardware/background processes were not strictly controlled.
- Observed: sensitivity report uses exhaustive Smith-Waterman as an internal reference ranking, not biological ground truth.
- Observed: cached exhaustive runtime may be reused instead of regenerated.
- Observed: benchmark docs mention current k-mer method as scan-based; indexed search has since appeared in source and tests, but the benchmark reports may not reflect indexed multi-search.
- Explicit limitation from reports: one query and score-derived tiers limit generalization.

## Tests that cover this behavior

- `tests/benchmarks/test_benchmark_alignment.py` checks benchmark smoke behavior and expected residue totals.
- Search benchmark scripts do not appear to have dedicated tests in the current tree.

## Questions or risks to verify manually

- See `architecture/open_questions.md` for benchmark report freshness relative to indexed search and cached-result provenance.
