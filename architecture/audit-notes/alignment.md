# Alignment

## Purpose

Document current pairwise alignment behavior for Needleman-Wunsch global alignment, Smith-Waterman local alignment, shared scoring helpers, substitution-matrix use, traceback, and returned result structures.

## Location

- `bioseq/alignment/needleman_wunsch.py`
- `bioseq/alignment/smith_waterman.py`
- `bioseq/alignment/scoring.py`
- `bioseq/alignment/substitution_matrices.py`
- `bioseq/alignment/alignment_stats.py`

## Public entry points

- `global_alignment(s1, s2, match=1, mismatch=-1, gap_penalty=-2, matrix=None, return_all=False, structured=True)`
- `local_alignment(s1, s2, match=2, mismatch=-1, gap_penalty=-2, matrix=None, return_all=False, structured=True)`
- `get_best_scores(s1, s2, gap_penalty, matrix=None, match=2, mismatch=-1)`
- `score_pair(a, b, loaded_matrix=None, match_score=2, mismatch_score=-1)`
- `get_alignment_stats(algn_1, algn_2)` and lower-level stats helpers.

## Inputs

- `s1`, `s2`: sequence strings. Local alignment rejects empty strings; global alignment has no explicit empty-input guard in the current code.
- `match`, `mismatch`: numeric simple scoring values when no substitution matrix is used.
- `gap_penalty`: numeric penalty. Positive inputs are converted to negative values by `normalize_gap`.
- `matrix`: `None`, matrix name string such as `"BLOSUM62"`, or a loaded matrix-like object accepted by `score_pair`.
- `return_all`: `bool` controlling whether traceback explores all recorded best moves.
- `structured`: `bool` controlling dictionary output versus legacy `list[tuple[str, str]]`.

## Input normalization and validation

- `score_pair` uppercases residues before scoring.
- `get_scoring_matrix` loads a named matrix through Biopython or returns a supplied matrix object unchanged.
- `normalize_gap` returns `-abs(gap)`.
- `local_alignment` raises `ValueError` when either input sequence is empty.
- Alignment stats helpers raise `ValueError` when aligned strings have unequal lengths.
- There is no observed alphabet validation before matrix scoring; unsupported matrix symbols would be handled by the matrix object, not by repository code.

## Main algorithm or workflow

Needleman-Wunsch:

1. Normalize `gap_penalty`.
2. Load optional scoring matrix.
3. Initialize a full dynamic-programming grid with cumulative gap penalties in the first row and column.
4. Fill the grid from top-left to bottom-right using diagonal, left, and up candidates.
5. Build a movement matrix with all tied best moves recorded as lists containing `"diag"`, `"left"`, and/or `"up"`.
6. Trace back from bottom-right to top-left. When `return_all=False`, use the first recorded move at each step.
7. Return either structured metadata and per-alignment statistics or legacy alignment tuple list.

Smith-Waterman:

1. Reject empty input sequences.
2. Normalize `gap_penalty`.
3. Load optional scoring matrix.
4. Initialize zero-filled score grid and movement grid.
5. Fill grid with `max(diagonal, horizontal, vertical, 0)`.
6. Track all positive best-score positions.
7. When `return_all=False`, trace only the first best position.
8. Stop traceback when the score grid reaches zero.
9. Return structured metadata and per-alignment statistics, or legacy tuple list.

`get_best_scores` is a score-only Smith-Waterman helper. It fills a score grid and returns `(best_score, best_positions)` without reconstructing alignments.

## Data structures created or consumed

- Scoring grid: `list[list[int | float]]`, dimensions `(len(s1) + 1) x (len(s2) + 1)`.
- Global movement matrix: `list[list[list[str] | None]]`, with move labels `"diag"`, `"left"`, and `"up"`.
- Local movement matrix: `list[list[list[str | None]]]`, with move labels `"diagonal"`, `"left"`, `"vertical"`, and sometimes `None` for all-negative movement candidates.
- Alignment tuple in legacy output: `tuple[str, str]` containing aligned sequence 1 and aligned sequence 2.
- Structured global result dictionary:
  - `algorithm`: `"Needleman-Wunsch"`.
  - `mode`: `"global"`.
  - `sequence_1`, `sequence_2`: original inputs.
  - `score`: final grid score.
  - `scoring`: scoring metadata.
  - `num_alignments`: number of reconstructed alignments.
  - `alignments`: `list[dict]` of aligned strings plus stats.
- Structured local result dictionary:
  - `algorithm`: `"Smith-Waterman"`.
  - `mode`: `"local"`.
  - `sequence_1`, `sequence_2`: original inputs.
  - `score`: best local score.
  - `scoring`: scoring metadata.
  - `best_positions`: `list[tuple[int, int]]`.
  - `num_alignments`: number of reconstructed alignments.
  - `alignments`: `list[dict]` of aligned strings plus stats.
- Alignment stats dictionary:
  - `alignment_length`: `int`.
  - `matches`: `int`, non-gap matching columns.
  - `mismatches`: `int`, non-gap differing columns.
  - `gaps`: `int`, total gap characters across both aligned strings.
  - `gap_columns`: `int`, columns containing at least one gap.
  - `identity`: rounded fraction including gaps.
  - `identity_excluding_gaps`: rounded fraction after subtracting gap characters from length.
  - `similarity`: currently `None`.

## Return value

- `global_alignment` returns a structured dictionary by default, or `list[tuple[str, str]]` when `structured=False`.
- `local_alignment` returns a structured dictionary by default, or `list[tuple[str, str]]` when `structured=False`.
- `get_best_scores` returns `(best_score, best_positions)`.

## Side effects

- Alignment functions do not read or write files.
- `load_matrix` uses `functools.lru_cache`, so matrix loading has process-local cache state.

## Dependencies called

- `Bio.Align.substitution_matrices.load`
- `bioseq.alignment.scoring`
- `bioseq.alignment.alignment_stats`

## Assumptions and limitations

- Observed: scoring metadata reports `"gap_model": "linear"`.
- Observed: no affine gap state matrices or separate open/extend penalties are present.
- Observed: `similarity` is always `None` in alignment statistics.
- Observed: matrix-based scoring omits `match` and `mismatch` fields from returned scoring metadata.
- Inferred implementation detail: when `return_all=False`, traceback chooses the first move in each stored tie list; this affects which optimal alignment is returned.
- Explicit limitation from README and benchmark docs: no affine gap penalties.

## Tests that cover this behavior

- `tests/alignment/test_needleman_wunsch_structured.py`
- `tests/alignment/test_smith_waterman.py`
- `tests/alignment/test_alignment_stats.py`
- `tests/alignment/test_scoring.py`
- selected CLI alignment tests in `tests/test_cli.py`.

## Questions or risks to verify manually

- See `architecture/open_questions.md` for empty global-alignment input, unsupported symbols with substitution matrices, and traceback tie-order risks.
