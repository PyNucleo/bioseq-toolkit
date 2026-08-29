# Search pipeline

## Public contracts

```python
search(query, database, k=3, threshold=1, top_n_hits=10,
       refinement=False, match_score=1, mismatch_score=-1,
       gap_penalty=-2, matrix=None)

multi_search(query_fasta, database, k=3, threshold=1, top_n_hits=10,
             indexed=True, refinement=False, match_score=1,
             mismatch_score=-1, gap_penalty=-2, matrix=None)
```

`database` is required. Query FASTA input may be a string path, `pathlib.Path`,
or query records. Database
normalization accepts `SequenceDatabase`, `list[str]`, or a string FASTA path.

## Validation and ordering

Both APIs validate before mode selection: `k` and `threshold` must be
non-Boolean integers at least 1. Search has no general alphabet compatibility
validation and no stable explicit `top_n_hits` type/range-validation contract.

The candidate sequence is:

1. generate unique uppercase k-mers;
2. apply the absolute threshold;
3. apply fixed `shared_kmers >= 0.3 * best_count` filtering;
4. rank by descending `shared_kmers`, then ascending `id`;
5. select `top_n_hits`;
6. optionally refine that selected set and re-rank by descending `sw_score`.

Refinement neither restores removed candidates nor performs another top-N
truncation. Regular mode runs `search()` per query. Indexed mode, the default,
builds one presence-based database-wide index and applies equivalent per-query
filtering/ranking. The index stores sequence membership, not seed positions.

## Scoring and results

With `refinement=False`, scoring options are inactive. With `matrix=None`,
Smith-Waterman uses match, mismatch, and a linear gap penalty. Matrix mode
replaces match/mismatch residue scores while keeping the linear gap penalty;
simple values remain accepted but inactive.

Candidate hits contain `id`, `sequence`, and `shared_kmers`; refined hits add
`sw_score` and `best_positions`. Each multi-query entry is:

```python
{"query_id": ..., "query_sequence": ..., "query_hits": [...]}
```

Python `best_positions` uses tuples; JSON serialization produces arrays.
Query order is preserved. Tests cover validation, ratio-filter equivalence,
deterministic ties, default/indexed behavior, both-mode refinement/scoring,
duplicate-ID rejection, and CLI JSON output.

## Selective property-based equivalence coverage

Hypothesis is a development/test dependency, used selectively for invariants
that benefit from generated combinations and shrinking; it is not Bioseq's
default testing mechanism. The first property lives in
`tests/pipelines/test_indexed_multi_search_equivalence.py` and checks the
candidate-search invariant that regular and indexed **unrefined** multi-search
produce exactly the same observable result. Indexed mode is a faster candidate
retrieval path for the same caller-visible operation, so a different answer is
not a successful optimization.

The property directly generates one to ten query records and one to ten raw
database sequences. Query IDs are arbitrary generated text of length 3–10 and
unique within each generated query list; query and database sequences use
canonical DNA characters `ATGCatgc` with lengths 4–10. Duplicate query
sequences and duplicate database sequence contents are allowed. It varies
`k` from 1–8, `threshold` from 1–4, and `top_n_hits` from 1–4, then calls
`multi_search(..., indexed=False, refinement=False)` and
`multi_search(..., indexed=True, refinement=False)` and asserts full Python
structure equality. Since the candidate stage has deterministic descending
shared-k-mer ranking with ascending-ID tie-breaking, this equality includes
observable query and hit ordering as well as returned values; it is not merely
set or hit-ID equivalence. The settings are `max_examples=150` and
`deadline=None`: the former bounds exploration rather than proving the
contract exhaustively, while the latter prevents runtime variation from being
reported as a correctness failure.

This intentionally bounded property does not cover proteins, ambiguous DNA,
malformed records or databases, duplicate query IDs, or refinement. Disabling
refinement isolates regular/indexed candidate retrieval; adding
Smith–Waterman scoring, refined ordering, and post-refinement behavior would
make any divergence harder to attribute and could mask a candidate-stage
error. Refined regular/indexed equivalence is a possible future invariant.

Deterministic pytest coverage remains the source for known examples, exact
boundaries, and named regressions, including ratio-filter parity,
sort-before-`top_n_hits`, deterministic ID tie-breaking, and a concrete
multi-query equality case. The Hypothesis property complements those tests
with generated interaction coverage and shrinking; it does not replace them.

This is software-correctness evidence over the configured generated domain,
not biological validation, search sensitivity or homology evidence,
Smith–Waterman/refinement validation, runtime evidence, an indexed speedup
claim, production-readiness evidence, BLAST equivalence, or exhaustive proof.
It is the correctness layer ahead of a separate future regular-versus-indexed
performance benchmark, whose datasets, parameters, timing method, environment,
repetitions, and limitations must be evaluated independently. Hypothesis test
runtime is not benchmark evidence.

## Limits

This is an educational filter/refine pipeline, not BLAST. It has no seed
extension, E-values, bit scores, affine gaps, significance calibration, or
general biological validation. No current benchmark establishes that indexed
mode is faster than regular scanning.
