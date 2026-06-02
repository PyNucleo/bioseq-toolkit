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
- Optional Smith-Waterman refinement of k-mer search hits
- Basic sequence database normalization
- Smith-Waterman runtime benchmarking on protein FASTA datasets
- k-mer-only search benchmarking across multiple `k` and threshold settings
- k-mer + Smith-Waterman refinement benchmarking
- Case-insensitive k-mer generation
- Benchmark driver scripts for alignment and search benchmarks
- Unit tests for core utilities, alignment, search, refinement, database normalization, FASTA loading, and benchmark smoke checks
- Minimal command-line interface for search, local alignment, and global alignment
- Editable installation support through `pyproject.toml`


The project has real structure, tests, benchmark reports, and a coherent search-pipeline direction. However, it is still early-stage and biologically incomplete.

Current important limitations include:

- No affine gap penalties yet
- No E-values or bit scores yet
- No indexed k-mer search yet
- No seed-extension step yet
- Command-line interface is still minimal and early-stage
- No completed biological case study yet
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

Runtime dependency:

- Biopython

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

The exact ordering depends on the candidate scores and filtering behavior.

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
        "sw_score": 10,
        "best_positions": [(5, 5)]
    }
]
```

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

Substitution-matrix support currently uses Biopython matrix loading. This is an initial implementation and should not yet be treated as a fully optimized scoring backend.

---

## Repository Structure

```text
bioseq-toolkit/
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
│   │   ├── kmer_search.py
│   │   ├── refinement.py
│   │   └── similarity_search.py
│   │
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
│   ├── BENCHMARKS.md
│   ├── benchmark_alignment.py
│   ├── benchmark_search.py
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
│   ├── test_cli.py
│   ├── test_fasta_io.py
│   ├── test_main.py
│   └── test_sequence_utils.py
│
├── pyproject.toml
├── pytest.ini
├── requirements.txt
└── README.md
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

Current limitation:

- UniProt-style and generic FASTA headers are supported.
- Other specialized formats such as RefSeq, PDB, and ASTRAL-specific headers may be added later when needed.

---

### Database Normalization

The search pipeline accepts:

- an existing `SequenceDatabase`
- a list of sequence strings
- a FASTA file path

The input is normalized into a consistent sequence-record format before search.

Example list-input normalization:

```python
[
    {"id": "id1", "sequence": "ATGCGT"},
    {"id": "id2", "sequence": "ATGCGA"}
]
```

For FASTA input, parsed FASTA metadata is preserved where available.

Current limitation:

- This is a lightweight wrapper, not a real database engine.
- It is intended to support consistent input handling for the current educational pipeline.

---

### Pairwise Alignment

The project currently includes two classic dynamic programming alignment algorithms:

- **Needleman-Wunsch** for global alignment
- **Smith-Waterman** for local alignment

Both support structured output.

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
- Matrix scoring exists initially, but scoring performance and broader validation still need improvement.

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

Current limitation:

- Matrix loading should be optimized further so matrices are loaded once and reused during dynamic programming.
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

The k-mer search step splits the query sequence and database sequences into words of length `k`.

It then:

1. generates unique query k-mers,
2. generates k-mers for each database sequence,
3. counts shared k-mers,
4. filters candidates by threshold,
5. applies a relative score filter,
6. returns candidate hits.

K-mer generation normalizes sequences to uppercase, so matching is case-insensitive for lowercase or mixed-case FASTA inputs.

This is the first step toward BLAST-like search behavior: use a fast word-based filter before doing more expensive alignment work.

Current limitation:

- The current implementation scans the database directly.
- It does not yet build an inverted k-mer index.
- It does not yet track seed positions.
- It does not yet perform seed extension.
- It does not currently validate whether the query and database are the same biological sequence type, such as DNA-vs-protein.

---

### Smith-Waterman Refinement

The search pipeline can optionally refine k-mer hits using Smith-Waterman local alignment.

Pipeline:

1. Normalize database input.
2. Run k-mer search.
3. Rank candidates by shared k-mers.
4. Keep the top candidate hits.
5. Optionally re-rank them using Smith-Waterman score.

This makes the project more than just a pairwise alignment implementation. It becomes a basic search pipeline demonstrating a seed-filter-refine idea.

Current limitation:

- Refinement currently scores top candidates, but does not yet return full local alignment objects inside each search hit.
- Search accuracy/sensitivity has not yet been fully evaluated against exhaustive Smith-Waterman rankings.

---

## Benchmarks

Benchmarking results are documented in:

```text
benchmarks/BENCHMARKS.md
```

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

These benchmarks demonstrate runtime improvement for the current query, datasets, implementation, and parameters.

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
- Current k-mer search scans database sequences directly rather than using an index
- Sensitivity and false-negative behavior have not yet been fully evaluated
- Current benchmarks emphasize speed more than biological accuracy

Future benchmarks should evaluate both:

1. runtime, and
2. candidate recovery compared with exhaustive Smith-Waterman search.

---

## Running Tests

Run all tests:

```bash
pytest
```

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
bioseq search
bioseq align-local
bioseq align-global
```

Run k-mer search:

```bash
bioseq search -q ATGCG -d data/benchmark_sequences/astral_10.fasta -k 3 -t 1
```

Run k-mer search with Smith-Waterman refinement:

```bash
bioseq search -q ATGCG -d data/benchmark_sequences/astral_10.fasta -k 3 -t 1 -r
```

Run Smith-Waterman local alignment:

```bash
bioseq align-local -s1 HEART -s2 HPEART --matrix BLOSUM62 -g -4
```

Run Needleman-Wunsch global alignment:

```bash
bioseq align-global -s1 ATGCG -s2 ATCGA -m 1 --mismatch -1 -g -2
```

The CLI prints structured JSON output.

The module-style form also remains supported:

```bash
python -m bioseq.cli --help
python -m bioseq.cli search -q ATGCG -d data/benchmark_sequences/astral_10.fasta
```

---

## Current Strengths

The strongest parts of the project are:

- Clear movement from isolated functions toward a search pipeline
- Basic package organization
- Structured output for both global and local alignment
- Structured FASTA parsing for UniProt-style and generic headers
- Case-insensitive k-mer generation for lowercase or mixed-case FASTA inputs
- Tests for alignment, search, refinement, database normalization, FASTA parsing, CLI behavior, and benchmark behavior
- Runtime comparison between exact dynamic programming and heuristic search
- Benchmark driver scripts for repeatable benchmark runs
- Minimal command-line interface for search, local alignment, and global alignment
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
- No indexed k-mer search yet
- Current k-mer search scans database sequences directly
- No seed-extension step yet
- Command-line interface is still minimal and early-stage
- No biological case study has been completed yet
- FASTA parsing currently supports structured UniProt-style and generic header metadata, but broader FASTA format support is limited
- K-mer search is case-insensitive, but it does not yet validate whether the query and database are the same biological sequence type, such as DNA-vs-protein
- Matrix scoring exists initially, but needs cleaner optimization and broader validation
- Not intended for production biological analysis

These limitations are intentional development targets, not hidden assumptions.

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
- Load substitution matrices once instead of repeatedly during scoring
- Add direct tests for known BLOSUM62 pair scores
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

### 4. Add Sensitivity / Recovery Evaluation

Future benchmarks should evaluate not only speed, but also whether heuristic search recovers the same important candidates as exhaustive Smith-Waterman.

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

### 6. Add Indexed Search

Build an inverted k-mer index so database k-mers do not need to be regenerated for every query.

Example idea:

```python
{
    "ATG": ["seq1", "seq5", "seq9"],
    "TGC": ["seq1", "seq2"]
}
```

This would make repeated searches more scalable.

---

### 7. Add Affine Gap Penalties

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

### 8. Add a Biological Case Study

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

### 9. Improve CLI and Packaging

A minimal command-line interface now exists for:

```bash
bioseq search
bioseq align-local
bioseq align-global
```

The module-style form also remains supported:

```bash
python -m bioseq.cli search
python -m bioseq.cli align-local
python -m bioseq.cli align-global
```

Current CLI output is structured JSON.

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
indexed search
        ↓
real protein-family case study
```

---

## License

No license has been specified yet.