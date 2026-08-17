# Translation

## Public contract

`process_fasta_sequences(file_path)` reads records through the shared strict
FASTA parser and returns:

```python
{"accepted": [...], "rejected": [...], "summary": {...}}
```

Accepted records contain `id`, one-based `record_position`, uppercase
`sequence`, `length`, `gc_content`, `transcribed_strand`, and a plain-string
`amino_acid_chain`. Rejected records contain `id`, `record_position`, uppercase
`sequence`, `reason_code`, `reason`, one-based `invalid_positions`, and a sorted
unique `invalid_symbols` list.

`summary` derives `total_records`, `accepted_records`, and `rejected_records`
from the two outcome lists. Record positions keep duplicate FASTA IDs
attributable to their source records.

## Workflow and boundaries

Each parsed sequence is uppercased. Symbols outside `A/T/G/C` are explicitly
rejected rather than silently skipped. Accepted DNA is measured, transcribed
through `mrna_template`, translated, and converted to a JSON-safe string.
Malformed FASTA structure propagates as `ValueError`; an empty FASTA produces
empty lists and zero summary counts.

The low-level `translate_sequence()` may return a Biopython `Seq`, while the
pipeline returns strings. Mapping utilities may raise `KeyError` when called
directly with unsupported symbols, but the pipeline diagnoses those symbols
before transcription. This is a limited DNA-symbol contract, not general
biological quality validation or ambiguous-base support.
