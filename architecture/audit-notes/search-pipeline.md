# Search pipeline

## Purpose

Document the current single-query and multi-query search behavior. The observed implementation is a k-mer candidate search pipeline with optional Smith-Waterman refinement for single-query search and regular multi-search. Indexed multi-search exists, but indexed refinement is explicitly rejected with `NotImplementedError`.

## Location

- `bioseq/pipelines/search_pipeline.py`
- `bioseq/search/kmer_search.py`
- `bioseq/search/kmer_index.py`
- `bioseq/search/similarity_search.py`
- `bioseq/search/refinement.py`

## Public entry points

- `search(query, database=None, k=3, threshold=1, top_n_hits=10, refinement=False)`
- `multi_search(query_fasta, database=None, k=3, threshold=1, top_n_hits=10, indexed=True, refinement=False, sort_results=True)`
- Supporting public helpers used directly by tests: `kmer_search`, `generate_kmers`, `validate_kmer_params`, `filter_by_relative_score`, `multi_query_indexed_search`, `refine_hits`, `rank_by_shared_kmers`.

## Inputs

- `query`: `str` sequence for `search`.
- `query_fasta`: either a FASTA path (`str` or `pathlib.Path`) or a list of query record dictionaries for `multi_search`.
- `database`: `SequenceDatabase`, `list[str]`, or `str` FASTA path accepted through `normalize_database`.
- `k`: `int` k-mer length.
- `threshold`: `int` minimum shared-k-mer count.
- `top_n_hits`: `int` maximum number of ranked hits returned per query.
- `refinement`: `bool`, controls Smith-Waterman re-ranking where supported.

## Input normalization and validation

- `normalize_database` converts `list[str]` into `SequenceDatabase([{"id": "id1", "sequence": seq1}, ...])`.
- `normalize_database` returns an existing `SequenceDatabase` unchanged.
- `normalize_database` treats `str` as a FASTA path and parses it through `database.load_database.create_database`.
- `multi_search` reads query records with `read_fasta_records` when `query_fasta` is `str` or `Path`; otherwise it uses the object directly.
- `validate_kmer_params` raises `TypeError` when `k` or `threshold` is not an integer, raises `ValueError` when `k <= 0` or `threshold < 0`, and returns `False` when `k > len(seq)`.
- `generate_kmers` uppercases sequences before creating unique k-mers.
- No biological alphabet validation is performed in search code; protein-like and DNA-like strings are both accepted by current tests.

## Main algorithm or workflow

Single-query `search`:

1. Normalize `database` to `SequenceDatabase`.
2. Run `kmer_search(query, db, k, threshold)`.
3. `kmer_search` builds a `set[str]` of query k-mers, scans every database record, counts shared unique k-mers per record, applies `threshold`, tracks the maximum shared-k-mer count, then applies `filter_by_relative_score(max_kmers, results, ratio=0.3)`.
4. `rank_by_shared_kmers` converts each candidate into a new hit dictionary and sorts by `(-shared_kmers, id)`.
5. The ranked list is sliced to `top_n_hits`.
6. If `refinement=True`, `refine_hits` adds `sw_score` and `best_positions` using `get_best_scores`, then sorts by descending `sw_score`.

Regular `multi_search(indexed=False)`:

1. Normalize/read query records.
2. For each query record, call `search(query["sequence"], ...)`.
3. Return one result dictionary per query with `query_id`, `query_sequence`, and `query_hits`.

Indexed `multi_search(indexed=True)`:

1. Reject `refinement=True` with `NotImplementedError`.
2. Normalize database and build a database-wide inverted k-mer index once.
3. For each query, generate unique query k-mers and count matching indexed database IDs.
4. Apply threshold inside indexed hit generation.
5. Apply the same 0.3 relative-score filter in `run_indexed_multi_search`.
6. Sort by shared-kmer count and ID through `rank_by_shared_kmers`.
7. Slice to `top_n_hits` per query.

## Data structures created or consumed

- `SequenceDatabase.sequences`: list stored on the object. In search tests this is `list[dict[str, str]]`, where each record has:
  - `id`: sequence identifier.
  - `sequence`: biological sequence string.
- Query records consumed by `multi_search`: `list[dict[str, str]]` with:
  - `id`: query identifier.
  - `sequence`: query sequence string.
- K-mer set: `set[str]` of uppercase unique k-mers, produced by `generate_kmers` and consumed by `get_shared_kmers` or `get_word_occurenes`.
- Candidate hit dictionaries returned by `kmer_search` and indexed search:
  - `id`: database record ID.
  - `sequence`: database sequence string.
  - `shared_kmers`: `int` count of unique query k-mers also present in the database sequence.
- Indexed search structure returned by `index_database_words`:
  - `kmer_index`: `dict[str, set[str]]`, mapping each k-mer to database record IDs that contain it.
  - `sequence_lookup`: `dict[str, str]`, mapping database record ID to sequence.
- Multi-search result dictionaries:
  - `query_id`: query record ID.
  - `query_sequence`: query sequence string.
  - `query_hits`: `list[dict]` of hit dictionaries.
- Refined hit dictionaries add:
  - `sw_score`: best Smith-Waterman local alignment score.
  - `best_positions`: list of matrix positions returned by `get_best_scores`; tuples in Python API, JSON arrays through CLI serialization.

## Return value

- `search` returns `list[dict]` sorted by shared k-mers and ID when `refinement=False`, or by descending `sw_score` when `refinement=True`.
- `multi_search` returns `list[dict]`, one result per query, preserving query order from the input records.
- Empty candidate sets return `[]`.

## Side effects

- Search itself has no file writes.
- When inputs are FASTA paths, search reads files through `read_fasta_records`.
- CLI use prints JSON to stdout.

## Dependencies called

- `database.database_utils.normalize_database`
- `bioseq.fasta_io.read_fasta_records`
- `bioseq.search.kmer_search.kmer_search`
- `bioseq.search.kmer_search.filter_by_relative_score`
- `bioseq.search.kmer_index.multi_query_indexed_search`
- `bioseq.search.similarity_search.rank_by_shared_kmers`
- `bioseq.search.refinement.refine_hits`
- `bioseq.alignment.smith_waterman.get_best_scores`

## Assumptions and limitations

- Observed: k-mer search uses unique k-mers, not positional seeds and not repeated k-mer frequencies.
- Observed: threshold filtering and relative-score filtering both affect returned candidates.
- Observed: `rank_by_shared_kmers` sorts equal scores by ID; tests describe this as deterministic tie-breaking, not biological evidence.
- Observed: indexed search is presence-based and database-wide, but not position-based.
- Observed: indexed multi-search refinement is not implemented.
- Inferred implementation detail: `sort_results` is accepted by `run_indexed_multi_search`, `run_regular_multi_search`, and `multi_search`, but current code does not use it to change behavior.
- Explicit limitation from benchmark docs: benchmark speedups are measurement artifacts for current query/data/parameters, not proof of broad biological performance.

## Tests that cover this behavior

- `tests/pipelines/test_search_pipeline.py`
- `tests/pipelines/test_search_many.py`
- `tests/pipelines/test_indexed_multi_search_equivalence.py`
- `tests/search/test_kmer_search.py`
- `tests/search/test_refinement.py`
- `tests/database/test_fasta_database_loading.py`
- CLI coverage for selected `search` commands in `tests/test_cli.py`.

## Questions or risks to verify manually

- See `architecture/open_questions.md` for `sort_results`, CLI `multi-search`, and indexed-refinement questions.
