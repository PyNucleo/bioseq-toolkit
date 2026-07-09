# Open questions

These items are unresolved because the current repository evidence is incomplete, inconsistent, untested, or environment-blocked. They should be verified manually before treating the behavior as a stable public guarantee.

## 1. CLI multi-search argument names

- Relevant file and symbol: `bioseq/cli.py`, `main`, `multi-search` branch.
- Why uncertain: the parser defines `--query-sequences`, `--database`, `--kmer-size`, `--threshold`, and `--regular`, but the execution branch references `args.top_n_hits`, `args.regular_kmer_search`, and possibly `args.query_sequences`. `top_n_hits` and `regular_kmer_search` are not defined by the observed parser.
- Manual check: run `python -m bioseq.cli multi-search ...` with a query FASTA and database FASTA, then add or inspect tests for this command.

## 2. `sort_results` parameter behavior

- Relevant file and symbol: `bioseq/pipelines/search_pipeline.py`, `multi_search`, `run_indexed_multi_search`, `run_regular_multi_search`.
- Why uncertain: `sort_results` is accepted and passed through but not used to alter sorting in the observed code.
- Manual check: decide whether `sort_results=False` is intended to preserve raw order, then add tests or remove/document the parameter later.

## 3. Indexed multi-search refinement

- Relevant file and symbol: `bioseq/pipelines/search_pipeline.py`, `run_indexed_multi_search`.
- Why uncertain: current code explicitly raises `NotImplementedError` for `refinement=True`, and tests confirm the error. There is no implementation-level design for indexed refinement.
- Manual check: if indexed refinement is desired, define whether refinement should occur before or after `top_n_hits` and relative-score filtering.

## 4. FASTA writer simple-header accession selection

- Relevant file and symbol: `bioseq/fasta_io.py`, `write_fasta_records`.
- Why uncertain: the code checks `"accession" in records`, but `records` is used as an iterable of record dictionaries. Current tests pass because `id` equals `accession` in the simple-header case, not because accession selection is proven.
- Manual check: add a test where `record["id"] != record["accession"]` and `full_header=False`.

## 5. UniProt FASTA response validation

- Relevant file and symbol: `bioseq/fasta_io.py`, `fetch_uniprot_sequences`; `case_studies/uniprot_fetch_demo/README.md`.
- Why uncertain: case-study documentation says fetched responses are treated as valid only when they contain FASTA text beginning with `>`, but the code only checks that HTTP 200 response text is non-empty before parsing.
- Manual check: mock `requests.get` with a non-empty non-FASTA body and inspect the returned record behavior.

## 6. `requests` dependency classification

- Relevant file and symbol: `pyproject.toml`, `requirements.txt`, `bioseq/fasta_io.py`.
- Why uncertain: `bioseq/fasta_io.py` imports `requests` at module import time, but `pyproject.toml` lists `requests` only under optional `dev` and `fetch` dependencies, not runtime dependencies.
- Manual check: install only the base project dependencies and import `bioseq.fasta_io`.

## 7. Global alignment empty-input behavior

- Relevant file and symbol: `bioseq/alignment/needleman_wunsch.py`, `global_alignment`.
- Why uncertain: local alignment explicitly rejects empty inputs; global alignment does not. Tests cover local empty rejection but not global empty behavior.
- Manual check: run `global_alignment("", "ATGC")`, `global_alignment("ATGC", "")`, and define intended behavior.

## 8. Substitution-matrix unsupported symbols

- Relevant file and symbol: `bioseq/alignment/scoring.py`, `score_pair`; `bioseq/alignment/substitution_matrices.py`, `load_matrix`.
- Why uncertain: score inputs are uppercased but not validated before indexing a matrix. Behavior for unsupported symbols depends on the Biopython matrix object.
- Manual check: call local/global alignment with `matrix="BLOSUM62"` and residues outside the matrix alphabet.

## 9. Traceback tie ordering as output contract

- Relevant files and symbols: `bioseq/alignment/needleman_wunsch.py`, `trace`; `bioseq/alignment/smith_waterman.py`, `local_trace`.
- Why uncertain: when `return_all=False`, current code follows the first stored movement. Tests verify selected examples but do not state whether first-move tie ordering is intended public behavior.
- Manual check: identify tied alignments and decide whether deterministic first alignment should be documented as stable or treated as implementation detail.

## 10. `SequenceDatabase` mutation and invariants

- Relevant file and symbol: `database/sequence_database.py`, `SequenceDatabase.get_sequences`.
- Why uncertain: `get_sequences` returns the internal list directly, and constructor input is stored without validation.
- Manual check: decide whether direct mutation is intended, and whether records must always contain `id` and `sequence`.

## 11. `normalize_database` path handling

- Relevant file and symbol: `database/database_utils.py`, `normalize_database`.
- Why uncertain: database paths are accepted only as `str`, while query FASTA paths in `multi_search` accept `str` and `Path`.
- Manual check: pass a `pathlib.Path` database path to `search` or `multi_search` and decide whether support should be added or documented as unsupported.

## 12. Translation result types

- Relevant file and symbol: `bioseq/translation.py`, `translate_sequence`; `bioseq/pipelines/translation_pipeline.py`, `process_fasta_sequences`.
- Why uncertain: `translate_sequence` returns a Biopython `Seq` object, while tests compare pipeline output to string literals. This works under current Biopython equality behavior but may confuse downstream callers.
- Manual check: decide whether public docs should promise `Seq` objects or strings.

## 13. Translation output key naming

- Relevant file and symbol: `bioseq/pipelines/translation_pipeline.py`, `process_fasta_sequences`.
- Why uncertain: the returned key `base_counts` contains GC percentages from `gc_content`, not base-count dictionaries from `base_number`.
- Manual check: decide whether this is intended naming, a docs issue, or a future code change.

## 14. Silent skipping in translation pipeline

- Relevant file and symbol: `bioseq/pipelines/translation_pipeline.py`, `process_fasta_sequences`.
- Why uncertain: empty and invalid-DNA records are skipped without reporting which records were skipped.
- Manual check: decide whether silent skipping is intended for educational use or whether skipped-record reporting is needed later.

## 15. Benchmark report freshness relative to indexed search

- Relevant files: `benchmarks/BENCHMARKS.md`, `benchmarks/SEARCH_SENSITIVITY.md`, `bioseq/search/kmer_index.py`, `tests/pipelines/test_indexed_multi_search_equivalence.py`.
- Why uncertain: benchmark reports describe the current k-mer method as scan-based and mention indexed search as future or next work, while source and tests now include indexed multi-search.
- Manual check: decide whether benchmark reports should be regenerated or annotated after indexed search work.

## 16. Cached benchmark result provenance

- Relevant files: `benchmarks/benchmark_search_sensitivity.py`, `benchmarks/cached_exhaustive_sw_results/astral_10000_BLOSUM62_gap-4.csv`, `benchmarks/SEARCH_SENSITIVITY.md`.
- Why uncertain: the sensitivity benchmark can reuse cached exhaustive Smith-Waterman runtime and scores. The cache is useful, but the docs do not make the regeneration command and exact provenance fully machine-verifiable.
- Manual check: record the exact command and environment used to produce the cached CSV if this benchmark will be cited externally.

## 17. Pytest temp-directory permissions in current environment

- Relevant command: `python -m pytest`; retried as `python -m pytest --basetemp .pytest-tmp -p no:cacheprovider`.
- Why uncertain: both attempts were blocked by `PermissionError` during pytest temporary-directory setup or cleanup. The observed failures were environmental setup errors, not project assertion failures, but a full clean test result was not obtained in this session.
- Manual check: fix local Windows temp/workspace permissions and rerun the full suite.
