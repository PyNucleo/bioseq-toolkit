# FASTA I/O

## Shared structural parser

`read_fasta_records(file_path)` opens local paths with a context manager and
passes lines to `_parse_fasta_records`. Nonempty downloaded UniProt response
bodies use the same parser. It:

- ignores blank lines;
- strips only CR/LF endings;
- rejects sequence before the first header;
- rejects empty IDs and records with empty sequences;
- rejects whitespace inside nonblank sequence lines;
- concatenates multiline sequence text;
- provides line-aware `ValueError` messages where applicable; and
- returns `[]` for empty or blank input.

String and path-like local paths accepted by `open()` work. Header routing
supports UniProt `sp|`/`tr|` and generic metadata. The parser validates FASTA
structure, not sequence alphabets.

## UniProt fetch

`fetch_uniprot_sequences(accession_file, strict=False)` reads one accession per
nonempty line and returns a flat `{"records": [...], "failed": [...]}` result.
Expected HTTP failures, `requests` exceptions, and empty HTTP-200 bodies are
recorded and skipped in non-strict mode; strict mode raises on the first such
operational failure. A nonempty malformed body reaches the shared parser, whose
`ValueError` propagates even when `strict=False`.

## Writer

`write_fasta_records(records, output_path, full_header)` requires every record
to contain a sequence. Full mode requires a usable stored header. Short mode
chooses per record: prefer a nonempty accession, otherwise use `id`; strip any
leading `>` before writing one marker. Missing sequence/header/identifier fields
raise `ValueError` with the one-based record number. Output uses one sequence
line per record and does not perform alphabet validation or wrapping.
