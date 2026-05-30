# bioseq-toolkit

Educational bioinformatics toolkit implementing sequence alignment, FASTA handling, k-mer search, and BLAST-like search ideas from scratch.

This project is built as a learning-focused bioinformatics software project. It is not intended to replace production tools such as BLAST, EMBOSS, Biopython, or professional sequence-analysis pipelines. Instead, the goal is to implement core sequence-analysis algorithms directly, test them, benchmark them, and gradually build toward a biologically realistic sequence-search toolkit.

## Current Project Status

`bioseq-toolkit` currently supports:

- Basic DNA sequence utilities
- FASTA parsing
- DNA transcription helpers
- Translation using Biopython
- Needleman-Wunsch global alignment
- Smith-Waterman local alignment
- k-mer based sequence search
- Optional Smith-Waterman refinement of k-mer search hits
- Basic sequence database normalization
- Smith-Waterman runtime benchmarking on protein FASTA datasets
- Unit tests for core utilities, search, refinement, database normalization, alignment, and benchmark smoke tests

The project is currently best described as an **educational sequence-search prototype**. It has real structure, tests, and benchmarks, but it is still early-stage and biologically incomplete.

## Why This Project Exists

Sequence alignment and database similarity search are central ideas in bioinformatics. Tools such as BLAST are fast because they avoid running full dynamic programming against every possible sequence unless needed. This project is an attempt to build those ideas step by step:

1. Start with basic sequence manipulation.
2. Implement exact pairwise alignment algorithms.
3. Build a simple k-mer search filter.
4. Refine promising hits using Smith-Waterman alignment.
5. Benchmark exact search behavior.
6. Gradually improve biological realism and scalability.

The long-term goal is to turn this into a small but coherent educational toolkit for understanding how sequence search works internally.

## Installation

Clone the repository:

```bash
git clone https://github.com/PyNucleo/bioseq-toolkit.git
cd bioseq-toolkit

```

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the test suite:

```bash
pytest
```

## Quick Start

Example: run a simple k-mer search.

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

Example with Smith-Waterman refinement:

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
    refinement=True
)

print(results)
```

When refinement is enabled, candidate hits are first found using shared k-mers, then re-ranked using Smith-Waterman local alignment scores.

## Repository Structure

```text
bioseq-toolkit/
├── bioseq/
│   ├── alignment/
│   │   ├── needleman_wunsch.py
│   │   ├── smith_waterman.py
│   │   └── scoring.py
│   ├── pipelines/
│   │   ├── search_pipeline.py
│   │   └── translation_pipeline.py
│   ├── search/
│   │   ├── kmer_search.py
│   │   ├── refinement.py
│   │   └── similarity_search.py
│   ├── fasta_io.py
│   ├── sequence_utils.py
│   ├── translation.py
│   └── validators.py
│
├── database/
│   ├── sequence_database.py
│   ├── database_utils.py
│   └── load_database.py
│
├── benchmarks/
│   ├── benchmark_alignment.py
│   ├── benchmark_utils.py
│   └── BENCHMARKS.md
│
├── data/
│   └── benchmark_sequences/
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
│   └── search/
│
├── requirements.txt
└── README.md
```

## Main Components

### Sequence Utilities

The project includes basic DNA sequence utility functions such as:

- sequence length
- base counting
- GC content
- reverse complement
- transcription from DNA template strand
- transcription from DNA coding strand

These utilities are intentionally simple and are tested as foundational building blocks.

### FASTA Parsing

The FASTA reader loads sequence records from FASTA files and returns structured records containing:

- generated sequence ID
- original FASTA header
- sequence string

This allows later search and benchmark code to work with a consistent record format.

### Pairwise Alignment

The project currently includes two classic dynamic programming alignment algorithms:

- **Needleman-Wunsch** for global alignment
- **Smith-Waterman** for local alignment

These implementations are educational and currently use a simple scoring system rather than full biological substitution matrices.

### k-mer Search

The k-mer search step splits the query sequence and database sequences into words of length `k`. It then counts shared k-mers and returns candidate hits that pass a threshold.

This is the first step toward BLAST-like search behavior: use a fast word-based filter before doing more expensive alignment work.

### Smith-Waterman Refinement

The search pipeline can optionally refine k-mer hits using Smith-Waterman local alignment.

Pipeline:

1. Normalize database input.
2. Run k-mer search.
3. Rank candidates by shared k-mers.
4. Keep the top candidate hits.
5. Optionally re-rank them using Smith-Waterman score.

This makes the project more than just a pairwise alignment implementation; it becomes a basic search pipeline.

## Benchmarks

Benchmarking results are documented in:

```text
benchmarks/BENCHMARKS.md
```

Current benchmarks evaluate Smith-Waterman runtime scaling on ASTRAL/SCOPe protein sequence datasets.

The benchmark report currently compares:

- score-only Smith-Waterman computation
- exhaustive alignment reconstruction
- runtime scaling across increasing dataset sizes
- total residue count as an explanation for runtime growth

The benchmarked datasets include chunks of:

- 10 sequences
- 100 sequences
- 1000 sequences
- 10000 sequences

The benchmark report also documents current limitations, including:

- linear gap model only
- simple scoring system
- no affine gap penalties yet
- exact Smith-Waterman benchmark only
- no current k-mer search benchmark yet

## Running Tests

Run all tests:

```bash
pytest
```

Run only search pipeline tests:

```bash
pytest tests/pipelines/
```

Run only benchmark tests:

```bash
pytest tests/benchmarks/
```

The tests currently cover:

- DNA utility functions
- database normalization
- Needleman-Wunsch alignment
- Smith-Waterman alignment
- k-mer search
- search refinement
- full search pipeline behavior
- benchmark smoke tests
- benchmark residue-count checks

## Examples

A simple search example is provided in:

```text
examples/search_demo.py
```

Run it with:

```bash
python examples/search_demo.py
```

## Current Strengths

The strongest parts of the project are:

- Clear movement from isolated algorithms toward a search pipeline
- Basic package organization
- Tests for search, refinement, database normalization, and benchmark behavior
- Benchmark documentation instead of only toy examples
- Honest educational scope
- Beginning of reproducibility through dataset chunks and benchmark reports

The project is especially useful for learning how sequence database search can be built from smaller algorithmic pieces.

## Current Limitations

This project is still early-stage. Important limitations include:

- No BLOSUM or PAM substitution matrix support yet
- No affine gap penalties yet
- No statistical significance estimates such as E-values or bit scores
- No indexed k-mer search yet
- Current k-mer search scans database sequences directly
- Alignment result objects are still minimal
- No stable command-line interface yet
- No biological case study has been completed yet
- Not intended for production biological analysis

These limitations are intentional development targets, not hidden assumptions.

## Roadmap

Planned development stages:

### 1. Improve Documentation and Reproducibility

- Expand README documentation
- Keep benchmark results synchronized with benchmark scripts
- Add clearer instructions for running examples and benchmarks

### 2. Improve Alignment Output

- Return structured alignment result objects
- Include alignment score
- Include aligned sequences
- Include identity percentage
- Include gap counts
- Include local alignment start/end positions where possible

### 3. Add Biologically Realistic Scoring

- Add BLOSUM62 support for protein alignment
- Add support for additional scoring matrices later
- Add affine gap penalties with separate gap opening and gap extension costs

### 4. Benchmark Search Strategies

Compare:

- exhaustive Smith-Waterman search
- k-mer search only
- k-mer search with Smith-Waterman refinement

Measure:

- runtime
- dataset size
- total residue count
- top-hit recovery
- speed/sensitivity tradeoffs

### 5. Add Indexed Search

Build an inverted k-mer index so that database k-mers do not need to be regenerated for every query.

Example idea:

```python
{
    "ATG": ["seq1", "seq5", "seq9"],
    "TGC": ["seq1", "seq2"]
}
```

This would make repeated searches more scalable.

### 6. Add a Biological Case Study

Use the toolkit on a real protein-family dataset.

Possible case study direction:

- collect related protein sequences
- search against a local FASTA database
- refine hits with Smith-Waterman
- compare conserved regions
- discuss biological interpretation and limitations

## Educational Scope

This repository is meant to show the internal logic behind common bioinformatics sequence-analysis tasks.

It is not currently designed for clinical, diagnostic, or production research use.

For serious biological analysis, established tools such as BLAST, HMMER, EMBOSS, MAFFT, MUSCLE, Clustal Omega, and Biopython should be used. This project is mainly for learning, experimentation, and building algorithmic understanding.

## Dependencies

Current dependencies are listed in:

```text
requirements.txt
```

Main dependencies:

- Biopython
- pytest
- pandas

## License

No license has been specified yet.