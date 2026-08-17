# Benchmarks

Benchmark artifacts are historical measurements for one implementation, query,
dataset family, parameter set, and environment. They are not general biological
performance evidence.

## Measurement families

The hard-coded protein query is 151 residues.

- `benchmark_alignment.py` measures full Smith-Waterman reconstruction and
  score-only `get_best_scores` with simple scoring `match=2`, `mismatch=-1`,
  `gap=-2`.
- `benchmark_search.py` measures scan-based `kmer_search`; its refined driver
  calls public `search()` with defaults `match=1`, `mismatch=-1`, `gap=-2`.
- `benchmark_search_sensitivity.py` measures scan-based k-mer retrieval and
  refinement against an internal exhaustive-SW ranking using BLOSUM62 and
  `gap=-4`.
- `benchmarks/main.py` only prints static historical values; it does not execute
  benchmark functions.

The sensitivity cache can reuse exhaustive scores and a stored runtime. Reference
tiers are derived from Smith-Waterman scores, not curated homology, family, or
function labels. Exhaustive Smith-Waterman is an internal comparator, not
biological ground truth.

## Current interpretation

Current toolkit source supports indexed multi-search, but these historical
measurements are scan-based and contain no regular-versus-indexed benchmark.
They do not establish indexed performance superiority. The one-query ASTRAL
results must be phrased as applying only to that query, dataset, scoring, and
recorded run.

Open work includes reproducible commands/environment capture, cache provenance,
regular-versus-indexed timing, harder and biologically curated cases, and affine
gap experiments. Recorded timing values are retained rather than silently
rewritten when implementation capabilities change.
