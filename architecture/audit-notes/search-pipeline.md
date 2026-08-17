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

## Limits

This is an educational filter/refine pipeline, not BLAST. It has no seed
extension, E-values, bit scores, affine gaps, significance calibration, or
general biological validation. No current benchmark establishes that indexed
mode is faster than regular scanning.
