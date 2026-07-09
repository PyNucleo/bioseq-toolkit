# FASTA I/O

## Purpose

Document current FASTA parsing, header interpretation, UniProt fetching, and FASTA writing behavior.

## Location

- `bioseq/fasta_io.py`
- `case_studies/uniprot_fetch_demo/README.md`
- `case_studies/uniprot_fetch_demo/accessions.txt`

## Public entry points

- `read_fasta_records(file)`
- `parse_fasta_header(header)`
- `parse_uniprot_fasta(header)`
- `parse_generic_fasta(header)`
- `read_fasta_sequences_only(FILE)`
- `parse_header_sequence_from_string(text)`
- `fetch_uniprot_sequences(accession_file, strict=False)`
- `write_fasta_records(records, output_path, full_header)`

## Inputs

- FASTA file path accepted by `open(file, "r")`.
- FASTA header string, normally beginning with `>`.
- UniProt accession file path: text file with one accession per non-empty line.
- `records` for writing: iterable of dictionaries with at least `sequence` and either `header`, `accession`, or `id` depending on `full_header`.
- `output_path`: path accepted by `open(output_path, "w", newline="")`.

## Input normalization and validation

- `read_fasta_records` strips trailing whitespace per line with `rstrip` and concatenates non-header lines as sequence text.
- `parse_fasta_header` strips the leading `>` for routing only.
- Headers beginning with `sp|` or `tr|` after stripping `>` are parsed as UniProt-style headers.
- Other headers are parsed as generic FASTA headers.
- `parse_uniprot_fasta` falls back to generic parsing if the split on `|` does not produce three parts.
- `fetch_uniprot_sequences` skips blank accession lines.
- `fetch_uniprot_sequences` treats HTTP 200 with empty response text as a failure; it does not currently check in code that response text begins with `>`.
- `write_fasta_records` does not validate line wrapping, sequence alphabet, or required keys before writing.

## Main algorithm or workflow

FASTA reading:

1. Open file.
2. Iterate lines.
3. On a header line, finalize the prior record if present.
4. Parse the new header into metadata and initialize `sequence`.
5. Concatenate sequence lines until the next header.
6. Append the final record.

UniProt header parsing:

1. Strip `>`.
2. Split into database prefix, accession, and remainder.
3. Use first remainder token as `entry_name`.
4. Collect description tokens until a token starting with `OS=`.

Generic header parsing:

1. Strip `>`.
2. Split on whitespace.
3. Use the first token as `id`.
4. Join remaining tokens into `description`.

UniProt fetching:

1. Read accession file.
2. For each non-empty accession, request `https://www.uniprot.org/uniprotkb/{accession}.fasta`.
3. For HTTP 200 with non-empty text, parse the first line as header and the remaining lines as sequence.
4. For non-200 responses, either raise in strict mode or append a failed-accession record.
5. Return both successful records and failures.

FASTA writing:

1. If `full_header=True`, write `record["header"]` exactly, then `record["sequence"]`.
2. If `full_header=False`, prepend `>` and write either a selected simple identifier or sequence ID, then sequence.

## Data structures created or consumed

- Parsed UniProt record dictionary:
  - `id`: accession string.
  - `db`: UniProt database prefix such as `"sp"` or `"tr"`.
  - `accession`: accession string.
  - `entry_name`: UniProt entry name, or empty string if absent.
  - `description`: description text before `OS=`.
  - `header`: original header string.
  - `sequence`: sequence string added by `read_fasta_records` or `fetch_uniprot_sequences`.
- Parsed generic record dictionary:
  - `id`: first token after `>`, or empty string for empty header.
  - `db`: `None`.
  - `accession`: `None`.
  - `entry_name`: `None`.
  - `description`: remaining header text.
  - `header`: original header string.
  - `sequence`: sequence string added by `read_fasta_records`.
- `fetch_uniprot_sequences` return dictionary:
  - `records`: `list[dict]` of parsed FASTA records.
  - `failed`: `list[dict]` with `accession`, `status_code`, and `reason`.
- `parse_header_sequence_from_string` return dictionary:
  - `header`: first line of stripped text.
  - `sequence`: all remaining lines concatenated.

## Return value

- Readers return parsed record lists or sequence-only lists.
- Header parsers return metadata dictionaries.
- Fetch returns a dictionary with successful and failed records.
- Writer returns `None`.

## Side effects

- `read_fasta_records` opens and reads local files.
- `fetch_uniprot_sequences` performs network I/O with `requests.get`.
- `write_fasta_records` writes local files.

## Dependencies called

- `requests.get`
- UniProt FASTA endpoint.

## Assumptions and limitations

- Observed: parser does not close files through a context manager in `read_fasta_records`.
- Observed: parser concatenates sequence lines without whitespace cleanup beyond `rstrip`.
- Observed: no sequence alphabet validation occurs in FASTA I/O.
- Observed: writer emits one sequence line per record and does not wrap FASTA output.
- Inferred implementation detail: `write_fasta_records(full_header=False)` checks `"accession" in records`, where `records` is expected to be a list; this likely does not detect per-record accession keys.
- Existing case-study notes describe stricter FASTA validity expectations than the current code fully enforces.

## Tests that cover this behavior

- `tests/test_fasta_io.py`
- `tests/test_fasta_writer.py`
- `tests/database/test_fasta_database_loading.py`
- selected CLI search and fetch command structure in `tests/test_cli.py`, though remote fetching itself is not covered by tests observed here.

## Questions or risks to verify manually

- See `architecture/open_questions.md` for simple-header accession selection, UniProt response validation, and fetch dependency configuration.
