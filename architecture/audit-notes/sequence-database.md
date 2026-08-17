# Sequence database

## Normalization boundary

`normalize_database(db)` accepts an existing `SequenceDatabase`, `list[str]`,
or a string FASTA path. Lists become records with generated IDs `id1`, `id2`,
and so on; FASTA metadata is preserved for string paths. `pathlib.Path` is not
accepted as a database object.

After every supported input form is normalized, duplicate IDs are rejected.
The `ValueError` lists every duplicated ID and all of its one-based record
positions. Identical sequences with different IDs are valid.

## Wrapper boundary

`SequenceDatabase` remains a lightweight mutable wrapper. Its constructor stores
the supplied object without enforcing record shape, duplicate-ID uniqueness, or
immutability, and `get_sequences()` returns that object directly. Therefore,
direct wrapper construction does not prove that public normalization invariants
were checked; search entry points call normalization before use.

Consumers expect records with `id` and `sequence`. There is no persistence or
database engine, and database `Path` support plus stronger mutation/shape
invariants remain open design questions.
