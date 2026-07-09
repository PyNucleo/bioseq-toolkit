# Sequence database

## Purpose

Document the minimal database abstraction currently used to pass sequence records across FASTA loading, search, indexed search, and benchmark code.

## Location

- `database/sequence_database.py`
- `database/database_utils.py`
- `database/load_database.py`

## Public entry points

- `SequenceDatabase(sequences=None)`
- `SequenceDatabase.get_sequences()`
- `normalize_database(db)`
- `create_database(db)`

## Inputs

- `SequenceDatabase` constructor accepts `None` or a supplied `sequences` object.
- `normalize_database` accepts:
  - existing `SequenceDatabase`;
  - `list`, interpreted as a list of sequence strings;
  - `str`, interpreted as a FASTA file path.
- `create_database` accepts a FASTA file path.

## Input normalization and validation

- `SequenceDatabase(None)` stores an empty list.
- `SequenceDatabase(sequences)` stores the object as provided; there is no type validation.
- `normalize_database(list)` creates records with generated IDs `id1`, `id2`, etc.
- `normalize_database(str)` calls `create_database`, which reads FASTA records and wraps them.
- Unsupported `normalize_database` input types raise bare `TypeError`.

## Main algorithm or workflow

1. Code that needs database records calls `normalize_database`.
2. Existing `SequenceDatabase` instances pass through unchanged.
3. Raw `list` inputs are converted into record dictionaries.
4. FASTA path strings are parsed into record dictionaries by `read_fasta_records`.
5. Search and index code call `db.get_sequences()` to retrieve the stored list.

## Data structures created or consumed

- `SequenceDatabase` class:
  - stored attribute `sequences`: usually `list[dict[str, str]]`.
  - known invariant from tests: when created from `list[str]`, duplicate sequences receive distinct generated IDs.
  - no enforced invariant for custom constructor input.
- Normalized record from `list[str]`:
  - `id`: generated string such as `"id1"`.
  - `sequence`: original sequence string.
- FASTA-derived record:
  - all keys produced by `read_fasta_records`, including at least `id`, `header`, and `sequence`.
- Consumers:
  - `kmer_search` expects records with `id` and `sequence`.
  - `index_database_words` expects records with `id` and `sequence`.
  - benchmark code expects FASTA-derived or normalized records depending on script.

## Return value

- `SequenceDatabase.get_sequences()` returns the stored `sequences` object directly.
- `normalize_database` and `create_database` return `SequenceDatabase`.

## Side effects

- `create_database` reads local FASTA files through `read_fasta_records`.
- `normalize_database` reads files when passed a `str`.
- No database persistence or mutation API is present beyond direct access to the stored list.

## Dependencies called

- `bioseq.fasta_io.read_fasta_records`
- `database.sequence_database.SequenceDatabase`

## Assumptions and limitations

- Observed: `SequenceDatabase` is a lightweight wrapper, not an indexed or validated database.
- Observed: `get_sequences` returns the internal object directly, so callers can mutate stored records.
- Observed: `normalize_database` treats all strings as file paths, not as single sequences.
- Observed: `Path` objects are not accepted by `normalize_database`; `multi_search` accepts `Path` only for query FASTA handling before database normalization.

## Tests that cover this behavior

- `tests/database/test_database_utils.py`
- `tests/database/test_fasta_database_loading.py`
- search and indexed-search tests that construct `SequenceDatabase` directly.

## Questions or risks to verify manually

- See `architecture/open_questions.md` for `Path` database support and mutation/invariant expectations.
