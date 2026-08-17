# Alignment and scoring

## Public entry points

- `global_alignment(s1, s2, match=1, mismatch=-1, gap_penalty=-2, matrix=None, return_all=False, structured=True)`
- `local_alignment(s1, s2, match=2, mismatch=-1, gap_penalty=-2, matrix=None, return_all=False, structured=True)`
- `get_best_scores(s1, s2, gap_penalty, matrix=None, match=2, mismatch=-1)`

Both algorithms normalize gaps with `-abs(gap)` and use a linear gap model.
Smith-Waterman rejects either empty input. The intended global double-empty
contract remains unresolved.

## Traceback

Needleman-Wunsch records tied moves in diagonal, left, then up order.
`return_all=True` retains every recorded optimal traceback leaf in deterministic
order; `return_all=False` follows the first movement. Completeness, ordering,
nested ties, counts, structured output, and legacy tuple output are regression-
tested. Deduplication of rendered alignments is not promised. Smith-Waterman
tracks positive best positions and tied local branches, but a universal public
ordering promise for every local tie configuration is not established.

## Structured schemas

Global results contain `algorithm`, `mode`, `sequence_1`, `sequence_2`, `score`,
`scoring`, `num_alignments`, and `alignments`; local results also contain
`best_positions`. Alignment entries contain `aligned_sequence_1`,
`aligned_sequence_2`, `alignment_length`, `matches`, `mismatches`, `gaps`,
`gap_columns`, `identity`, `identity_excluding_gaps`, and `similarity`.
`similarity` is currently `None`, not a computed biological metric.

## Matrix boundaries

Simple scoring accepts unrestricted characters. Matrix paths uppercase residues
and raise `ValueError` for covered unknown names or unsupported residues. Named
and already-loaded matrices may provide different diagnostic context; exact
message equality across paths is not guaranteed. Matrix mode omits inactive
match/mismatch metadata and keeps the linear gap penalty.

Neither alignment implementation provides affine gaps, E-values, bit scores,
or biological significance calibration.
