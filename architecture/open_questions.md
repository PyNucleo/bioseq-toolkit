# Open questions

Resolved implementation questions have been removed. These remaining items are
not stable public guarantees.

## 1. Global alignment empty-input contract

`global_alignment` produces gap-only results when exactly one sequence is empty,
but the double-empty structured path currently fails while computing alignment
identity. Tests do not establish the intended public contract. This branch does
not document the inconsistent double-empty behavior as supported.

## 2. Local-alignment tie ordering

Needleman-Wunsch `return_all=True` completeness and deterministic movement order
are regression-tested, as is deterministic first-path behavior. Smith-Waterman
records and explores ties, but identical public ordering guarantees for every
local-alignment tie configuration are not established.

## 3. `SequenceDatabase` mutation and record-shape invariants

The constructor stores supplied objects without validation and `get_sequences()`
returns the internal object directly. Public normalization checks duplicate IDs,
but direct wrapper construction does not enforce record shape or immutability.

## 4. Database `Path` support

`normalize_database` accepts string FASTA database paths but not `pathlib.Path`.
Multi-search accepts `Path` only for its query FASTA. Whether database `Path`
objects should be supported remains a design question.

## 5. Benchmark freshness and reproducibility

Recorded reports are historical and scan-based. Current indexed multi-search has
no regular-versus-indexed performance benchmark. Hardware, environment, exact
commands, and cached exhaustive-SW provenance are not fully machine-verifiable,
so results should remain workload- and revision-scoped until regenerated with a
stronger provenance record.
