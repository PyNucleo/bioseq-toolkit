# Translation

## Purpose

Document the current FASTA-to-translation pipeline and lower-level sequence utilities used by that pipeline.

## Location

- `bioseq/pipelines/translation_pipeline.py`
- `bioseq/translation.py`
- `bioseq/validators.py`
- `bioseq/sequence_utils.py`

## Public entry points

- `process_fasta_sequences(file_path)`
- `translate_sequence(seq, stop_at_stop=True, trim_partial=True)`
- `valid_dna(seq)`
- `total_length(sequence)`
- `base_number(sequence)`
- `gc_content(sequence)`
- `reverse_complement(sequence)`
- `mrna_template(sequence)`
- `mrna_coding(sequence)`

## Inputs

- `process_fasta_sequences`: FASTA file path.
- `translate_sequence`: nucleotide or RNA-like string accepted by `Bio.Seq.Seq`.
- Sequence utilities: sequence strings.

## Input normalization and validation

- `process_fasta_sequences` reads records with `read_fasta_records`.
- Each sequence is uppercased.
- Empty sequences are skipped.
- Sequences failing `valid_dna` are skipped.
- `valid_dna` uppercases and requires every character to be one of `A`, `T`, `G`, or `C`.
- `translate_sequence(trim_partial=True)` trims sequence length down to a multiple of 3 before Biopython translation.
- Sequence utilities uppercase before counting/transcription/complement operations.
- `reverse_complement`, `mrna_template`, and `mrna_coding` use dictionary lookups and will raise `KeyError` for unsupported bases.

## Main algorithm or workflow

1. Read FASTA records.
2. For each record, uppercase `record["sequence"]`.
3. Skip empty or non-DNA sequences.
4. Store accepted DNA sequence.
5. Transcribe accepted DNA through `mrna_template`, using template-strand complement rules.
6. Compute length and GC content.
7. Translate transcribed mRNA with `translate_sequence`, trimming partial codons by default and stopping at the first stop codon by default.
8. Return parallel lists inside one dictionary.

## Data structures created or consumed

- Input FASTA records from `read_fasta_records`: dictionaries with at least `sequence`.
- Final result dictionary:
  - `sequence`: `list[str]` of accepted uppercase DNA sequences.
  - `length`: `list[int]`, each item from `total_length`.
  - `base_counts`: `list[float]`, despite name, each item is GC percentage from `gc_content`.
  - `transcribed_strand`: `list[str]` of mRNA strings produced by `mrna_template`.
  - `amino_acid_chain`: `list[Bio.Seq.Seq]` in implementation; tests compare to strings and pass through Biopython equality behavior.
- `base_number` returns `collections.Counter`.

## Return value

- `process_fasta_sequences` returns the final dictionary of parallel lists.
- `translate_sequence` returns a Biopython `Seq` translation result.
- `valid_dna` returns `bool`.
- `gc_content` returns rounded percentage as `float`.

## Side effects

- `process_fasta_sequences` reads a local FASTA file.
- No writes or network calls occur in translation code.

## Dependencies called

- `bioseq.fasta_io.read_fasta_records`
- `Bio.Seq.Seq.translate`
- `collections.Counter`

## Assumptions and limitations

- Observed: translation pipeline only accepts DNA alphabet `A/T/G/C`; RNA input containing `U` is rejected at pipeline level.
- Observed: transcription uses template-strand mapping, not coding-strand mapping.
- Observed: key `base_counts` actually stores GC percentages, not base-count dictionaries.
- Observed: invalid and empty records are silently skipped.
- Inferred implementation detail: Biopython accepts the mRNA string produced by `mrna_template` and returns `Seq` objects whose equality can compare with strings in tests.

## Tests that cover this behavior

- `tests/pipelines/test_translation_pipeline.py`
- `tests/test_sequence_utils.py`

## Questions or risks to verify manually

- See `architecture/open_questions.md` for `base_counts` naming, silent skip behavior, and returned `Seq` versus string expectations.
