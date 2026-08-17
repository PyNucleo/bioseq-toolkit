# bioseq-toolkit

Educational bioinformatics toolkit implementing core sequence-analysis algorithms from scratch, including FASTA handling, DNA sequence utilities, pairwise alignment, k-mer search, Smith-Waterman refinement, and benchmarked BLAST-like search ideas.

This project is built as a learning-focused bioinformatics software project. It is **not** intended to replace production tools such as BLAST, EMBOSS, HMMER, MAFFT, MUSCLE, Clustal Omega, Biopython, or professional sequence-analysis pipelines.

The goal is to implement the internal logic behind common bioinformatics tasks directly, test the implementations, benchmark them, and gradually improve biological realism, reproducibility, and scalability.

---

## Current Project Status

`bioseq-toolkit` is currently best described as an **educational sequence-search prototype**.

It currently supports:

- Basic DNA sequence utilities
- FASTA parsing
- DNA transcription helpers
- Translation using Biopython
- Needleman-Wunsch global alignment
- Smith-Waterman local alignment
- Structured output for both global and local alignment
- Alignment statistics, including:
  - alignment length
  - matches
  - mismatches
  - gaps
  - gap columns
  - identity
  - identity excluding gaps
- Simple match/mismatch scoring
- Linear gap penalties
- Initial substitution-matrix support through Biopython, including BLOSUM62 usage in alignment tests
- k-mer based sequence search
- Regular and database-wide indexed multi-query search (indexed by default)
- Optional Smith-Waterman refinement in single-query and both multi-query modes
- Basic sequence database normalization
- Strict structural FASTA parsing shared by local reads and UniProt responses
- Accountable FASTA translation results with accepted, rejected, and summary sections
- Smith-Waterman runtime benchmarking on protein FASTA datasets
- k-mer-only search benchmarking across multiple `k` and threshold settings
- k-mer + Smith-Waterman refinement benchmarking
- Case-insensitive k-mer generation
- Benchmark driver scripts for alignment and search benchmarks
- Unit tests for core utilities, alignment, search, refinement, database normalization, FASTA loading, and benchmark smoke checks
- Five-command interface for search, multi-search, local/global alignment, and UniProt fetching
- Editable installation support through `pyproject.toml`


The project has real structure, tests, benchmark reports, and a coherent search-pipeline direction. However, it is still early-stage and biologically incomplete.

Current important limitations include:

- No affine gap penalties yet
- No E-values or bit scores yet
- No seed-extension step yet
- No biologically validated protein-family case study yet
- Not intended for clinical, diagnostic, or production biological analysis

---

## Why This Project Exists

Sequence alignment and database similarity search are central ideas in bioinformatics.

Tools such as BLAST are fast because they do not simply run full dynamic programming against every database sequence in the most expensive way possible. Instead, they use heuristic ideas such as word matching, candidate filtering, and refinement.

This project builds toward those ideas step by step:

1. Start with basic biological sequence manipulation.
2. Implement exact pairwise alignment algorithms.
3. Build a simple k-mer search filter.
4. Refine promising hits using Smith-Waterman local alignment.
5. Benchmark exact search, k-mer search, and k-mer + refinement.
6. Gradually improve scoring, biological realism, scalability, and reproducibility.

The long-term goal is to turn this into a small but coherent educational toolkit for understanding how sequence search works internally.

A good description of the project is:

> An educational BLAST-like sequence-search prototype demonstrating exact dynamic programming alignment, k-mer candidate filtering, Smith-Waterman refinement, and runtime tradeoffs.

A bad description would be:

> A replacement for BLAST.

This project is not that.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/PyNucleo/bioseq-toolkit.git
cd bioseq-toolkit
```

Install dependencies:

```bash
pip install -e ".[dev]"
```

Run the test suite:

```bash
pytest
```

Using `python -m pytest` is recommended because it runs tests through the active Python environment and is less likely to run into import-path issues than calling `pytest` directly in some setups.

---

## Dependencies

Runtime dependencies:

- Biopython
- requests

Development dependency:

- pytest

Dependencies are defined in `pyproject.toml`. The `requirements.txt` file is kept as a lightweight compatibility/development helper.

---

## Quick Start

### Example 1 — Simple k-mer search

```python
from bioseq.pipelines.search_pipeline import search

results = search(
    query="ATGCG",
    database=[
        "ATGCGT",
        "ATGCGA",
        "GGGGGG"
    ],
    k=3,
    threshold=1
)

print(results)
```

Expected style of output:

```python
[
    {"id": "id1", "sequence": "ATGCGT", "shared_kmers": 3},
    {"id": "id2", "sequence": "ATGCGA", "shared_kmers": 3}
]
```

Ties are deterministic: hits are ordered by descending `shared_kmers`, then
ascending record `id`.

---

### Example 2 — k-mer search with Smith-Waterman refinement

```python
from bioseq.pipelines.search_pipeline import search

results = search(
    query="ATGCG",
    database=[
        "ATGAAA",
        "ATGCGT",
        "GGGGGG"
    ],
    k=3,
    threshold=1,
    top_n_hits=3,
    refinement=True
)

print(results)
```

When refinement is enabled:

1. Candidate hits are first found using shared k-mers.
2. The top candidates are scored using Smith-Waterman local alignment.
3. Hits are re-ranked by Smith-Waterman score.

Expected style of refined output:

```python
[
    {
        "id": "id2",
        "sequence": "ATGCGT",
        "shared_kmers": 3,
        "sw_score": 5,
        "best_positions": [(5, 5)]
    },
    {
        "id": "id1",
        "sequence": "ATGAAA",
        "shared_kmers": 1,
        "sw_score": 3,
        "best_positions": [(3, 3)]
    }
]
```

Refinement operates only on the already filtered and selected `top_n_hits`
candidate set; it does not restore candidates removed by k-mer filtering.

---

### Example 3 — Smith-Waterman local alignment

```python
from bioseq.alignment.smith_waterman import local_alignment

result = local_alignment(
    "GATTACA",
    "TTAC",
    match=2,
    mismatch=-1,
    gap_penalty=-2,
    structured=True
)

print(result)
```

Structured output includes:

- algorithm name
- alignment mode
- input sequences
- score
- scoring metadata
- best matrix positions
- number of alignments
- aligned sequences
- alignment statistics

---

### Example 4 — Needleman-Wunsch global alignment

```python
from bioseq.alignment.needleman_wunsch import global_alignment

result = global_alignment(
    "ATGCG",
    "ATCGA",
    match=1,
    mismatch=-1,
    gap_penalty=-2,
    structured=True
)

print(result)
```

---

### Example 5 — Alignment with BLOSUM62

```python
from bioseq.alignment.smith_waterman import local_alignment

result = local_alignment(
    "HEART",
    "HPEART",
    gap_penalty=-4,
    matrix="BLOSUM62",
    structured=True
)

print(result)
```

Substitution-matrix support uses cached Biopython matrix loading. BLOSUM62 paths
and repository error boundaries are tested, but this is not a broadly validated
biological scoring backend.

---

## Repository Structure

```text
bioseq-toolkit/
├── architecture/
│   ├── README.md
│   ├── model.c4
│   ├── open_questions.md
│   ├── audit-notes/
│   └── views/
│       ├── system-overview.c4
│       ├── search-workflow.c4
│       ├── alignment-workflow.c4
│       ├── fasta-and-database-flow.c4
│       ├── translation-and-benchmarks.c4
│       ├── search-pipeline-detail.c4
│       ├── kmer-candidate-retrieval-detail.c4
│       ├── smith-waterman-detail.c4
│       ├── needleman-wunsch-detail.c4
│       ├── fasta-database-normalization-detail.c4
│       └── translation-detail.c4
│
├── bioseq/
│   ├── alignment/
│   │   ├── alignment_stats.py
│   │   ├── needleman_wunsch.py
│   │   ├── scoring.py
│   │   ├── smith_waterman.py
│   │   └── substitution_matrices.py
│   │
│   ├── pipelines/
│   │   ├── search_pipeline.py
│   │   └── translation_pipeline.py
│   │
│   ├── search/
│   │   ├── kmer_index.py
│   │   ├── kmer_search.py
│   │   ├── refinement.py
│   │   └── similarity_search.py
│   │
│   ├── __init__.py
│   ├── cli.py
│   ├── fasta_io.py
│   ├── main.py
│   ├── sequence_utils.py
│   ├── translation.py
│   └── validators.py
│
├── database/
│   ├── database_utils.py
│   ├── load_database.py
│   └── sequence_database.py
│
├── benchmarks/
│   ├── figures/
│   ├── BENCHMARKS.md
│   ├── SEARCH_SENSITIVITY.md
│   ├── benchmark_alignment.py
│   ├── benchmark_search.py
│   ├── benchmark_search_sensitivity.py
│   ├── benchmark_utils.py
│   ├── main.py
│   ├── run_alignment_benchmarks.py
│   └── run_search_benchmarks.py
│
├── data/
│   └── benchmark_sequences/
│       ├── astral_10.fasta
│       ├── astral_100.fasta
│       ├── astral_1000.fasta
│       └── astral_10000.fasta
│
├── case_studies/
│   └── uniprot_fetch_demo/
├── dataset_tools/
│   └── chunk_dataset.py
│
├── examples/
│   └── search_demo.py
│
├── tests/
│   ├── alignment/
│   ├── benchmarks/
│   ├── database/
│   ├── pipelines/
│   ├── search/
│   └── test_sequence_utils.py
│
├── .gitignore
├── README.md
├── pyproject.toml
├── pytest.ini
└── requirements.txt
```

---

## Main Components

### Sequence Utilities

The project includes basic DNA sequence utility functions such as:

- sequence length
- base counting
- GC content
- reverse complement
- transcription from DNA template strand
- transcription from DNA coding strand

These are intentionally simple foundational functions.

Current limitation:

- Direct utility calls are still mainly designed around uppercase DNA strings.

---

### FASTA Parsing

The FASTA reader loads sequence records from FASTA files and returns structured records.

For UniProt-style headers such as:

```text
>sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens
```

the parser extracts:

```python
{
    "id": "P69905",
    "db": "sp",
    "accession": "P69905",
    "entry_name": "HBA_HUMAN",
    "description": "Hemoglobin subunit alpha",
    "header": ">sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens",
    "sequence": "..."
}
```

For generic FASTA headers such as:

```text
>seq1 some description
```

the parser returns:

```python
{
    "id": "seq1",
    "db": None,
    "accession": None,
    "entry_name": None,
    "description": "some description",
    "header": ">seq1 some description",
    "sequence": "..."
}
```

Local reads and nonempty UniProt response bodies use the same structural
parser. It ignores blank lines, strips only CR/LF line endings, concatenates
multiline sequence data, and rejects sequence data before the first header,
empty IDs, whitespace inside nonblank sequence lines, and empty record
sequences. Empty or blank input returns an empty list. Errors include line
context where applicable, and local files are opened with a context manager.
Both strings and path-like objects accepted by `open()` work for local FASTA
reads. These checks validate FASTA structure, not biological alphabets.

`fetch_uniprot_sequences(accession_file, strict=False)` reads one accession per
nonempty line and returns a flat result:

```python
{"records": [...], "failed": [...]}
```

In non-strict mode, expected HTTP/request failures and empty HTTP-200 bodies
are recorded in `failed` and fetching continues. With `strict=True`, the first
such operational failure raises. A nonempty malformed response is passed to
the shared parser; its `ValueError` propagates even in non-strict mode.

`write_fasta_records()` requires a usable sequence on every record. Full-header
mode requires a usable stored `header`. Short-header mode chooses independently
for each record, preferring a nonempty `accession` and otherwise using `id`; a
leading `>` is removed before the writer constructs one short header marker.
Missing required fields raise `ValueError` with the record number.

UniProt-style and generic headers are supported. Other specialized formats may
be added only when a real dataset requires them.

### Translation pipeline

`process_fasta_sequences(path)` returns accountable, JSON-safe results:

```python
{
    "accepted": [
        {
            "id": "accepted",
            "record_position": 1,
            "sequence": "ATGGCC",
            "length": 6,
            "gc_content": 66.67,
            "transcribed_strand": "UACCGG",
            "amino_acid_chain": "YR"
        }
    ],
    "rejected": [
        {
            "id": "rejected",
            "record_position": 2,
            "sequence": "ATGN",
            "reason_code": "unsupported_dna_symbols",
            "reason": "Sequence contains unsupported DNA symbols.",
            "invalid_positions": [4],
            "invalid_symbols": ["N"]
        }
    ],
    "summary": {
        "total_records": 2,
        "accepted_records": 1,
        "rejected_records": 1
    }
}
```

`record_position` and `invalid_positions` are one-based. Accepted DNA is
uppercased; unsupported symbols are reported with deterministic unique symbol
lists instead of being silently skipped. Duplicate FASTA IDs remain
attributable through record positions, and malformed FASTA errors propagate.
The pipeline uses the `A/T/G/C` alphabet; this is not general biological quality
validation. The lower-level `translate_sequence()` may still return a Biopython
`Seq`, while the pipeline converts translations to plain strings.

---

### Database Normalization

The search pipeline accepts:

- an existing `SequenceDatabase`
- a list of sequence strings
- a FASTA file path supplied as `str`

The input is normalized into a consistent sequence-record format before search.

Example list-input normalization:

```python
[
    {"id": "id1", "sequence": "ATGCGT"},
    {"id": "id2", "sequence": "ATGCGA"}
]
```

For FASTA input, parsed FASTA metadata is preserved where available.

Public normalization rejects duplicate record IDs and reports all one-based
positions for every duplicated ID. Identical sequence strings remain valid
when their IDs differ, including generated records from `list[str]` input.
`pathlib.Path` is accepted for a multi-search query FASTA but is not currently
accepted as the database object; database `Path` support remains an open design
question.

Current limitation:

- `SequenceDatabase` is a lightweight mutable wrapper, not a database engine.
- Direct construction stores its input without enforcing record shape or the
  public normalization invariants; `get_sequences()` exposes the stored object.
- It is intended to support consistent input handling for the current educational pipeline.

---

### Pairwise Alignment

The project currently includes two classic dynamic programming alignment algorithms:

- **Needleman-Wunsch** for global alignment
- **Smith-Waterman** for local alignment

Both support structured output.

Needleman-Wunsch uses a linear gap model. With `return_all=True`, it retains
all tied optimal traceback leaves in deterministic recorded movement order;
with `return_all=False`, it follows the deterministic first movement. Tied
counts, ordering, legacy tuple output, and structured metadata are covered by
regression tests. Rendered alignments are not promised to be deduplicated. The
double-empty global-alignment case remains unresolved and is not documented as
supported. Smith-Waterman also uses a linear gap model and does not provide
calibrated biological significance.

Example structured fields:

```python
{
    "algorithm": "Smith-Waterman",
    "mode": "local",
    "sequence_1": "...",
    "sequence_2": "...",
    "score": 8,
    "scoring": {
        "match": 2,
        "mismatch": -1,
        "gap_penalty": -2,
        "matrix": None,
        "gap_model": "linear"
    },
    "best_positions": [(6, 4)],
    "num_alignments": 1,
    "alignments": [
        {
            "aligned_sequence_1": "TTAC",
            "aligned_sequence_2": "TTAC",
            "alignment_length": 4,
            "matches": 4,
            "mismatches": 0,
            "gaps": 0,
            "gap_columns": 0,
            "identity": 1.0,
            "identity_excluding_gaps": 1.0,
            "similarity": None
        }
    ]
}
```

Current limitations:

- Gap model is currently linear.
- Affine gap penalties are not implemented yet.
- Similarity is currently a placeholder field.
- Matrix scoring is tested for supported BLOSUM62 paths, while broader
  biological validation remains limited.

---

### Scoring

The project currently supports:

- simple match/mismatch scoring
- linear gap penalties
- initial substitution-matrix lookup through Biopython

Simple scoring example:

```python
match = 2
mismatch = -1
gap_penalty = -2
```

Matrix scoring example:

```python
matrix = "BLOSUM62"
gap_penalty = -4
```

Repository scoring paths raise `ValueError` for covered unknown matrix names
and unsupported residues. Named and already-loaded matrix paths can provide
different diagnostic context; callers should not depend on one identical error
string for every path. Matrix residue-pair scores replace simple match/mismatch
values, but the linear `gap_penalty` remains active.

Current limitation:

- Matrix loading is cached, but matrix-based benchmark behavior is not yet fully evaluated.
- PAM matrices and custom user-provided matrices are not yet documented as stable features.
- Affine gap penalties are not yet implemented.

---

### Alignment Statistics

Alignment results can include:

- alignment length
- matches
- mismatches
- total gaps
- gap columns
- identity
- identity excluding gaps
- similarity placeholder

This makes outputs easier to test, inspect, and compare.

---

### k-mer Search

The public APIs are:

```python
search(query, database, k=3, threshold=1, top_n_hits=10,
       refinement=False, match_score=1, mismatch_score=-1,
       gap_penalty=-2, matrix=None)

multi_search(query_fasta, database, k=3, threshold=1, top_n_hits=10,
             indexed=True, refinement=False, match_score=1,
             mismatch_score=-1, gap_penalty=-2, matrix=None)
```

`database` is required. `multi_search` supports the regular path with
`indexed=False` and the database-wide indexed path with `indexed=True`; indexed
mode is the default. Both modes support the same optional refinement controls.

Candidate processing is:

1. generate unique uppercase k-mers;
2. keep records meeting the absolute shared-k-mer `threshold`;
3. apply the fixed relative cutoff `shared_kmers >= 0.3 * best_count`;
4. sort by descending `shared_kmers`, then ascending `id`;
5. select `top_n_hits`;
6. optionally Smith-Waterman-refine that selected subset and re-rank it by
   descending `sw_score`.

K-mer generation normalizes sequences to uppercase, so matching is case-insensitive for lowercase or mixed-case FASTA inputs.

`k` and `threshold` must be integers (Boolean values do not count), and both
must be at least 1. There is no stable explicit public type/range-validation
contract for `top_n_hits`, and search does not validate general biological
alphabet compatibility. Candidate hits contain `id`, `sequence`, and
`shared_kmers`. Each multi-query result uses this outer schema:

```python
{"query_id": "...", "query_sequence": "...", "query_hits": [...]}
```

This is the first step toward BLAST-like search behavior: use a fast word-based filter before doing more expensive alignment work.

Current limitations:

- The regular path scans database sequences directly; the indexed multi-query
  path builds a presence-based database-wide index without seed positions.
- Current historical benchmarks measure the scan-based path and do not
  establish that indexed mode is faster.
- It does not yet track seed positions.
- It does not yet perform seed extension.
- It does not currently validate whether the query and database are the same biological sequence type, such as DNA-vs-protein.

---

### Smith-Waterman Refinement

The search pipeline can optionally refine k-mer hits using Smith-Waterman local alignment.

Pipeline:

1. Normalize database input.
2. Filter and rank candidates as described above.
3. Keep `top_n_hits` candidates.
4. Refine only that selected set and re-rank by descending Smith-Waterman score.

With `refinement=False`, all scoring options are inactive. With refinement and
`matrix=None`, Smith-Waterman uses `match_score`, `mismatch_score`, and the
linear `gap_penalty`. With a named matrix, residue-pair matrix scores replace
simple match/mismatch scoring; simple values remain accepted but inactive, and
the linear gap penalty remains active. Refined hits add `sw_score` and
`best_positions`. Refinement does not run a second top-N truncation or
reintroduce candidates removed earlier.

This makes the project more than just a pairwise alignment implementation. It becomes a basic search pipeline demonstrating a seed-filter-refine idea.

Current limitation:

- Refinement currently scores top candidates, but does not yet return full local alignment objects inside each search hit.
- The recorded sensitivity work is limited to one query, dataset family, and
  parameter set; it is not broad biological validation.

---

## Benchmarks

Benchmarking results are documented in:

```text
benchmarks/BENCHMARKS.md
benchmarks/SEARCH_SENSITIVITY.md
```

These are historical recorded measurements, not results regenerated for every
revision. The hard-coded query is 151 residues. Alignment runtime measurements
use simple scoring `match=2`, `mismatch=-1`, `gap=-2`; the historical refined
public-search driver uses the API defaults `match=1`, `mismatch=-1`, `gap=-2`;
and sensitivity measurements use BLOSUM62 with a linear gap penalty of `-4`.
The runtime reports use scan-based k-mer functions. Current indexed
multi-search support has no regular-versus-indexed performance benchmark.

Current benchmarks evaluate runtime scaling on ASTRAL/SCOPe protein sequence datasets.

The benchmark report currently compares:

- score-only Smith-Waterman computation
- Smith-Waterman alignment reconstruction
- k-mer-only search using multiple `k` and threshold settings
- k-mer + Smith-Waterman refinement
- runtime scaling across increasing dataset sizes
- total residue count as an explanation for runtime growth

The benchmarked datasets include chunks of:

- 10 sequences
- 100 sequences
- 1000 sequences
- 10000 sequences

Benchmark FASTA chunks are stored in:

```text
data/benchmark_sequences/
```

Current benchmark highlights:

- Exact Smith-Waterman score-only search on the 10000-sequence dataset took about 137.18 seconds.
- K-mer-only search on the same dataset took about 0.50–0.52 seconds depending on `k` and threshold.
- K-mer + Smith-Waterman refinement took about 0.67 seconds with `k=3`, `threshold=3`, and `top_n_hits=10`.
- K-mer + refinement showed a runtime speedup of about 203.73× over exact Smith-Waterman score-only search on the 10000-sequence dataset.

These historical measurements demonstrate runtime improvement for that query,
those datasets, that implementation revision, and those parameters.

They do **not** prove equal biological sensitivity or accuracy compared with exhaustive Smith-Waterman search.

---

## Benchmark Drivers

Benchmark driver scripts are provided in the `benchmarks/` directory.

Run the Smith-Waterman alignment benchmarks:

```bash
python -m benchmarks.run_alignment_benchmarks
```

Run the k-mer and k-mer + Smith-Waterman refinement benchmarks:

```bash
python -m benchmarks.run_search_benchmarks
```

These benchmark drivers should be run from the repository root so dataset paths resolve correctly.

Full benchmark runs may take a long time, especially on larger datasets and repeated iterations.

---

## Benchmark Limitations

Current benchmark limitations:

- Linear gap model only
- No affine gap penalties yet
- Simple scoring is still the main benchmark mode
- Initial substitution-matrix support exists, but matrix-based benchmark behavior is not yet fully evaluated
- Protein lengths vary across datasets
- Hardware and background processes were not strictly controlled
- Historical runtime and sensitivity measurements use scan-based candidate
  retrieval; indexed multi-search performance has not been compared.
- Sensitivity evidence is limited to one query and score-derived internal tiers.
- Current benchmarks emphasize speed more than biological accuracy

Future benchmarks should evaluate both:

1. runtime, and
2. candidate recovery compared with exhaustive Smith-Waterman search.

---

## Running Tests

Run all tests:

```bash
python -m pytest -q
```

At documentation synchronization base revision `3ce3f74`, the full suite passes:
`199 passed`. This count is revision-scoped rather than a permanent project
claim.

Run only alignment tests:

```bash
python -m pytest tests/alignment/
```

Run only search pipeline tests:

```bash
python -m pytest tests/pipelines/
```

Run only benchmark tests:

```bash
python -m pytest tests/benchmarks/
```

The tests currently cover:

- DNA utility functions
- lowercase DNA utility behavior
- translation pipeline behavior
- database normalization
- FASTA database loading
- FASTA header parsing
- Needleman-Wunsch alignment
- Needleman-Wunsch structured output
- Smith-Waterman alignment
- Smith-Waterman structured output
- BLOSUM62 usage in alignment tests
- k-mer search
- case-insensitive k-mer behavior
- search refinement
- full search pipeline behavior
- CLI subprocess behavior
- benchmark smoke tests
- benchmark residue-count checks

---

## Examples

A simple search example is provided in:

```text
examples/search_demo.py
```

Run it from the repository root:

```bash
python -m examples.search_demo
```

---

## Command-Line Interface

The project includes a minimal command-line interface.

After editable installation:

```bash
pip install -e ".[dev]"
```

Run:

```bash
bioseq --help
```

Available commands:

```bash
bioseq search --help
bioseq multi-search --help
bioseq align-local --help
bioseq align-global --help
bioseq fetch-uniprot --help
```

`search`, `multi-search`, `align-local`, and `align-global` write JSON to
stdout. `fetch-uniprot` writes a FASTA file and prints a human-readable summary;
it is the exception to the JSON-output pattern.

Run k-mer search with a tiny FASTA file:

```bash
printf ">seq1 best hit\nATGCGT\n>seq2 also good\nATGCGA\n>seq3 no hit\nGGGGGG\n" > tiny_search.fasta
bioseq search -q ATGCG -d tiny_search.fasta -k 3 -t 2
```

Run k-mer search with Smith-Waterman refinement:

```bash
bioseq search -q ATGCG -d tiny_search.fasta -k 3 -t 1 -r
```

Run Smith-Waterman local alignment:

```bash
bioseq align-local -s1 HEART -s2 HPEART --matrix BLOSUM62 -g -4
```

Run Needleman-Wunsch global alignment:

```bash
bioseq align-global -s1 ATGCG -s2 ATCGA -m 1 --mismatch -1 -g -2
```

Run indexed multi-search (the default) or regular multi-search:

```bash
bioseq multi-search -q queries.fasta -d tiny_search.fasta -k 3 -t 1
bioseq multi-search -q queries.fasta -d tiny_search.fasta -k 3 -t 1 --regular
```

Both modes support `--refine`, `--match`, `--mismatch`, `--gap-penalty`, and
`--matrix` and emit the same `query_hits` schema.

Fetch accessions from a plain-text file:

```bash
bioseq fetch-uniprot -f accessions.txt -o sequences.fasta
```

`--strict` stops on the first expected operational fetch failure;
`--full-header` preserves stored headers instead of short accession/ID headers;
and `--show-failed` prints the failures collected in non-strict mode.

Python and CLI defaults are not all identical: Python `search()` defaults to
`threshold=1`, while CLI `search` defaults to `3`; Python `local_alignment()`
defaults to `match=2`, while CLI `align-local` defaults to `1`; CLI
`multi-search` defaults to `threshold=1`.

The module-style form also remains supported:

```bash
python -m bioseq.cli --help
python -m bioseq.cli search -q ATGCG -d tiny_search.fasta -k 3 -t 2
python -m bioseq.cli multi-search -q queries.fasta -d tiny_search.fasta
```

---

## Current Strengths

The strongest parts of the project are:

- Clear movement from isolated functions toward a search pipeline
- Basic package organization
- Structured output for both global and local alignment
- Structured FASTA parsing for UniProt-style and generic headers
- Shared strict FASTA structure parsing for local and downloaded text
- Accountable translation outcomes with source-record positions
- Case-insensitive k-mer generation for lowercase or mixed-case FASTA inputs
- Tested regular/indexed multi-search equivalence and both-mode refinement
- Tests for alignment, search, refinement, database normalization, FASTA parsing, CLI behavior, and benchmark behavior
- Runtime comparison between exact dynamic programming and heuristic search
- Benchmark driver scripts for repeatable benchmark runs
- Five-command interface with machine-readable JSON for search/alignment commands
- Editable installation support through `pyproject.toml`
- `bioseq` console command entry point
- Benchmark documentation instead of only toy examples
- Initial substitution-matrix support
- Honest educational scope
- Beginning of reproducibility through dataset chunks, benchmark reports, CLI tests, pytest configuration, and package metadata

The project is especially useful for learning how sequence database search can be built from smaller algorithmic pieces.

---

## Current Limitations

This project is still early-stage.

Important limitations include:

- No affine gap penalties yet
- No statistical significance estimates such as E-values or bit scores
- Regular search scans database sequences; indexed multi-search is supported,
  but no current benchmark establishes its performance advantage
- No seed-extension step yet
- No biologically meaningful protein-family case study has been completed yet
- FASTA parsing currently supports structured UniProt-style and generic header metadata, but broader FASTA format support is limited
- K-mer search is case-insensitive, but it does not yet validate whether the query and database are the same biological sequence type, such as DNA-vs-protein
- Matrix scoring has tested BLOSUM62 paths but still needs broader validation
- Not intended for production biological analysis

These limitations are intentional development targets.

---

## Roadmap

Planned development stages:

### 1. Improve Benchmark Reproducibility

- Keep benchmark results synchronized with benchmark scripts
- Continue using explicit keyword arguments in benchmark function calls
- Add clearer instructions for running benchmark modes
- Save benchmark outputs in a consistent machine-readable format where useful
- Separate speed benchmarks from sensitivity/recovery benchmarks

---

### 2. Clean Up Substitution Matrix Support

- Keep BLOSUM62 support
- Keep substitution matrix loading cached and covered by tests
- Preserve direct tests for known BLOSUM62 pair scores and `ValueError` boundaries
- Benchmark simple scoring vs BLOSUM62 scoring
- Clearly document which matrices are stable and tested

---

### 3. Improve FASTA and Database Handling

Completed for current scope.

Current FASTA parsing supports:

- UniProt-style reviewed headers such as `sp|...`
- UniProt-style unreviewed headers such as `tr|...`
- Generic FASTA headers as a safe fallback
- Preserved original FASTA headers
- Structured metadata fields for IDs, accessions, entry names, descriptions, and sequences

Future work:

- Add more specialized parsers only when needed, such as RefSeq, PDB, or ASTRAL-specific formats
- Add more FASTA edge-case tests if new datasets require them

---

### 4. Extend Sensitivity / Recovery Evaluation

An initial one-query ASTRAL sensitivity report now compares heuristic results
with an internal exhaustive Smith-Waterman ranking. Future benchmarks should
extend this beyond the recorded query, data, parameters, and score-derived tiers.

Possible measurements:

- top-1 recovery
- top-5 recovery
- top-10 recovery
- overlap with exact Smith-Waterman rankings
- effect of `k` and threshold on missed candidates
- speed/sensitivity tradeoff

---

### 5. Improve Search Heuristics

- Add seed-position tracking
- Add seed-extension behavior around k-mer hits
- Compare seed-extension search against current k-mer-only search
- Continue measuring runtime and candidate recovery

---

### 6. Add Affine Gap Penalties

Current alignments use a linear gap model.

Future affine gap support should separate:

- gap opening penalty
- gap extension penalty

Example future scoring metadata:

```python
{
    "gap_model": "affine",
    "gap_open": -10,
    "gap_extend": -1
}
```

This would make alignments more biologically realistic.

---

### 7. Add a Biological Case Study

Use the toolkit on a real protein-family dataset.

Possible case study direction:

- collect related protein sequences from a reliable source
- build a local FASTA database
- search against the database
- refine hits with Smith-Waterman
- compare exact search vs k-mer search vs k-mer + refinement
- align top hits using BLOSUM62
- discuss conserved regions
- explain biological interpretation and limitations

Possible target:

- beta-lactamase / antibiotic resistance protein family

A biological case study is important because it would turn the project from a software-only exercise into a bioinformatics analysis project.

---

### 8. Improve CLI and Packaging

A minimal command-line interface now exists for:

```bash
bioseq search --help
bioseq multi-search --help
bioseq align-local --help
bioseq align-global --help
bioseq fetch-uniprot --help
```

The module-style form also remains supported:

```bash
python -m bioseq.cli search --help
python -m bioseq.cli align-local --help
python -m bioseq.cli align-global --help
python -m bioseq.cli fetch-uniprot --help
python -m bioseq.cli multi-search --help
```

Search and alignment commands emit JSON; `fetch-uniprot`
writes FASTA and prints a human-readable summary.

Current packaging support includes:

- `pyproject.toml`
- editable installation with `pip install -e ".[dev]"`
- `bioseq` console command entry point
- runtime dependency declaration
- development test dependency declaration

Future CLI/package improvements:

- Add optional output-file support
- Add translation command later
- Add benchmark command wrappers later
- Improve CLI examples and documentation
- Keep the CLI thin: parse user arguments, call existing tested functions, and print structured output

---

## Educational Scope

This repository is meant to show the internal logic behind common bioinformatics sequence-analysis tasks.

It is useful for:

- learning pairwise alignment
- understanding dynamic programming
- understanding local vs global alignment
- experimenting with k-mer based filtering
- seeing why heuristic search can be faster than exhaustive search
- practicing testing and benchmarking of bioinformatics code
- practicing how bioinformatics code can be exposed through a command-line interface

It is not currently designed for:

- clinical analysis
- diagnostic use
- production research pipelines
- replacing mature bioinformatics tools

For serious biological analysis, established tools such as BLAST, HMMER, EMBOSS, MAFFT, MUSCLE, Clustal Omega, and Biopython should be used.

---

## Project Philosophy

This project prioritizes:

- correctness over speed of development
- tested behavior over hidden assumptions
- honest limitations over exaggerated claims
- incremental biological realism
- reproducible benchmarking
- learning by implementing core ideas directly

The project should not pretend to be more mature than it is.

The intended direction is:

```text
basic sequence utilities
        ↓
pairwise alignment
        ↓
k-mer candidate search
        ↓
Smith-Waterman refinement
        ↓
benchmarking
        ↓
biologically realistic scoring
        ↓
regular-vs-indexed benchmarking
        ↓
real protein-family case study
```

---

## License

No license has been specified yet.
