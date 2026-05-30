[1mdiff --git a/README.md b/README.md[m
[1mindex 98d24b4..8a6464e 100644[m
[1m--- a/README.md[m
[1m+++ b/README.md[m
[36m@@ -1,16 +1,416 @@[m
 # bioseq-toolkit[m
[31m-Educational bioinformatics toolkit implementing sequence alignment and BLAST-like search algorithms from scratch.[m
[32m+[m
[32m+[m[32mEducational bioinformatics toolkit implementing sequence alignment, FASTA handling, k-mer search, and BLAST-like search ideas from scratch.[m
[32m+[m
[32m+[m[32mThis project is built as a learning-focused bioinformatics software project. It is not intended to replace production tools such as BLAST, EMBOSS, Biopython, or professional sequence-analysis pipelines. Instead, the goal is to implement core sequence-analysis algorithms directly, test them, benchmark them, and gradually build toward a biologically realistic sequence-search toolkit.[m
[32m+[m
[32m+[m[32m## Current Project Status[m
[32m+[m
[32m+[m[32m`bioseq-toolkit` currently supports:[m
[32m+[m
[32m+[m[32m- Basic DNA sequence utilities[m
[32m+[m[32m- FASTA parsing[m
[32m+[m[32m- DNA transcription helpers[m
[32m+[m[32m- Translation using Biopython[m
[32m+[m[32m- Needleman-Wunsch global alignment[m
[32m+[m[32m- Smith-Waterman local alignment[m
[32m+[m[32m- k-mer based sequence search[m
[32m+[m[32m- Optional Smith-Waterman refinement of k-mer search hits[m
[32m+[m[32m- Basic sequence database normalization[m
[32m+[m[32m- Smith-Waterman runtime benchmarking on protein FASTA datasets[m
[32m+[m[32m- Unit tests for core utilities, search, refinement, database normalization, alignment, and benchmark smoke tests[m
[32m+[m
[32m+[m[32mThe project is currently best described as an **educational sequence-search prototype**. It has real structure, tests, and benchmarks, but it is still early-stage and biologically incomplete.[m
[32m+[m
[32m+[m[32m## Why This Project Exists[m
[32m+[m
[32m+[m[32mSequence alignment and database similarity search are central ideas in bioinformatics. Tools such as BLAST are fast because they avoid running full dynamic programming against every possible sequence unless needed. This project is an attempt to build those ideas step by step:[m
[32m+[m
[32m+[m[32m1. Start with basic sequence manipulation.[m
[32m+[m[32m2. Implement exact pairwise alignment algorithms.[m
[32m+[m[32m3. Build a simple k-mer search filter.[m
[32m+[m[32m4. Refine promising hits using Smith-Waterman alignment.[m
[32m+[m[32m5. Benchmark exact search behavior.[m
[32m+[m[32m6. Gradually improve biological realism and scalability.[m
[32m+[m
[32m+[m[32mThe long-term goal is to turn this into a small but coherent educational toolkit for understanding how sequence search works internally.[m
 [m
 ## Installation[m
 [m
[31m-Clone repository:[m
[32m+[m[32mClone the repository:[m
 [m
 ```bash[m
 git clone https://github.com/PyNucleo/bioseq-toolkit.git[m
 cd bioseq-toolkit[m
 [m
[31m-Run:[m
[32m+[m[32m```[m
[32m+[m
[32m+[m[32mInstall dependencies:[m
 [m
 ```bash[m
 pip install -r requirements.txt[m
[31m-pytest[m
\ No newline at end of file[m
[32m+[m[32m```[m
[32m+[m
[32m+[m[32mRun the test suite:[m
[32m+[m
[32m+[m[32m```bash[m
[32m+[m[32mpytest[m
[32m+[m[32m```[m
[32m+[m
[32m+[m[32m## Quick Start[m
[32m+[m
[32m+[m[32mExample: run a simple k-mer search.[m
[32m+[m
[32m+[m[32m```python[m
[32m+[m[32mfrom bioseq.pipelines.search_pipeline import search[m
[32m+[m
[32m+[m[32mresults = search([m
[32m+[m[32m    query="ATGCG",[m
[32m+[m[32m    database=[[m
[32m+[m[32m        "ATGCGT",[m
[32m+[m[32m        "ATGCGA",[m
[32m+[m[32m        "GGGGGG"[m
[32m+[m[32m    ],[m
[32m+[m[32m    k=3,[m
[32m+[m[32m    threshold=1[m
[32m+[m[32m)[m
[32m+[m
[32m+[m[32mprint(results)[m
[32m+[m[32m```[m
[32m+[m
[32m+[m[32mExpected style of output:[m
[32m+[m
[32m+[m[32m```python[m
[32m+[m[32m[[m
[32m+[m[32m    {"id": "id1", "sequence": "ATGCGT", "shared_kmers": 3},[m
[32m+[m[32m    {"id": "id2", "sequence": "ATGCGA", "shared_kmers": 3}[m
[32m+[m[32m][m
[32m+[m[32m```[m
[32m+[m
[32m+[m[32mExample with Smith-Waterman refinement:[m
[32m+[m
[32m+[m[32m```python[m
[32m+[m[32mfrom bioseq.pipelines.search_pipeline import search[m
[32m+[m
[32m+[m[32mresults = search([m
[32m+[m[32m    query="ATGCG",[m
[32m+[m[32m    database=[[m
[32m+[m[32m        "ATGAAA",[m
[32m+[m[32m        "ATGCGT",[m
[32m+[m[32m        "GGGGGG"[m
[32m+[m[32m    ],[m
[32m+[m[32m    k=3,[m
[32m+[m[32m    threshold=1,[m
[32m+[m[32m    refinement=True[m
[32m+[m[32m)[m
[32m+[m
[32m+[m[32mprint(results)[m
[32m+[m[32m```[m
[32m+[m
[32m+[m[32mWhen refinement is enabled, candidate hits are first found using shared k-mers, then re-ranked using Smith-Waterman local alignment scores.[m
[32m+[m
[32m+[m[32m## Repository Structure[m
[32m+[m
[32m+[m[32m```text[m
[32m+[m[32mbioseq-toolkit/[m
[32m+[m[32m├── bioseq/[m
[32m+[m[32m│   ├── alignment/[m
[32m+[m[32m│   │   ├── needleman_wunsch.py[m
[32m+[m[32m│   │   ├── smith_waterman.py[m
[32m+[m[32m│   │   └── scoring.py[m
[32m+[m[32m│   ├── pipelines/[m
[32m+[m[32m│   │   ├── search_pipeline.py[m
[32m+[m[32m│   │   └── translation_pipeline.py[m
[32m+[m[32m│   ├── search/[m
[32m+[m[32m│   │   ├── kmer_search.py[m
[32m+[m[32m│   │   ├── refinement.py[m
[32m+[m[32m│   │   └── similarity_search.py[m
[32m+[m[32m│   ├── fasta_io.py[m
[32m+[m[32m│   ├── sequence_utils.py[m
[32m+[m[32m│   ├── translation.py[m
[32m+[m[32m│   └── validators.py[m
[32m+[m[32m│[m
[32m+[m[32m├── database/[m
[32m+[m[32m│   ├── sequence_database.py[m
[32m+[m[32m│   ├── database_utils.py[m
[32m+[m[32m│   └── load_database.py[m
[32m+[m[32m│[m
[32m+[m[32m├── benchmarks/[m
[32m+[m[32m│   ├── benchmark_alignment.py[m
[32m+[m[32m│   ├── benchmark_utils.py[m
[32m+[m[32m│   └── BENCHMARKS.md[m
[32m+[m[32m│[m
[32m+[m[32m├── data/[m
[32m+[m[32m│   └── benchmark_sequences/[m
[32m+[m[32m│[m
[32m+[m[32m├── dataset_tools/[m
[32m+[m[32m│   └── chunk_dataset.py[m
[32m+[m[32m│[m
[32m+[m[32m├── examples/[m
[32m+[m[32m│   └── search_demo.py[m
[32m+[m[32m│[m
[32m+[m[32m├── tests/[m
[32m+[m[32m│   ├── alignment/[m
[32m+[m[32m│   ├── benchmarks/[m
[32m+[m[32m│   ├── database/[m
[32m+[m[32m│   ├── pipelines/[m
[32m+[m[32m│   └── search/[m
[32m+[m[32m│[m
[32m+[m[32m├── requirements.txt[m
[32m+[m[32m└── README.md[m
[32m+[m[32m```[m
[32m+[m
[32m+[m[32m## Main Components[m
[32m+[m
[32m+[m[32m### Sequence Utilities[m
[32m+[m
[32m+[m[32mThe project includes basic DNA sequence utility functions such as:[m
[32m+[m
[32m+[m[32m- sequence length[m
[32m+[m[32m- base counting[m
[32m+[m[32m- GC content[m
[32m+[m[32m- reverse complement[m
[32m+[m[32m- transcription from DNA template strand[m
[32m+[m[32m- transcription from DNA coding strand[m
[32m+[m
[32m+[m[32mThese utilities are intentionally simple and are tested as foundational building blocks.[m
[32m+[m
[32m+[m[32m### FASTA Parsing[m
[32m+[m
[32m+[m[32mThe FASTA reader loads sequence records from FASTA files and returns structured records containing:[m
[32m+[m
[32m+[m[32m- generated sequence ID[m
[32m+[m[32m- original FASTA header[m
[32m+[m[32m- sequence string[m
[32m+[m
[32m+[m[32mThis allows later search and benchmark code to work with a consistent record format.[m
[32m+[m
[32m+[m[32m### Pairwise Alignment[m
[32m+[m
[32m+[m[32mThe project currently includes two classic dynamic programming alignment algorithms:[m
[32m+[m
[32m+[m[32m- **Needleman-Wunsch** for global alignment[m
[32m+[m[32m- **Smith-Waterman** for local alignment[m
[32m+[m
[32m+[m[32mThese implementations are educational and currently use a simple scoring system rather than full biological substitution matrices.[m
[32m+[m
[32m+[m[32m### k-mer Search[m
[32m+[m
[32m+[m[32mThe k-mer search step splits the query sequence and database sequences into words of length `k`. It then counts shared k-mers and returns candidate hits that pass a threshold.[m
[32m+[m
[32m+[m[32mThis is the first step toward BLAST-like search behavior: use a fast word-based filter before doing more expensive alignment work.[m
[32m+[m
[32m+[m[32m### Smith-Waterman Refinement[m
[32m+[m
[32m+[m[32mThe search pipeline can optionally refine k-mer hits using Smith-Waterman local alignment.[m
[32m+[m
[32m+[m[32mPipeline:[m
[32m+[m
[32m+[m[32m1. Normalize database input.[m
[32m+[m[32m2. Run k-mer search.[m
[32m+[m[32m3. Rank candidates by shared k-mers.[m
[32m+[m[32m4. Keep the top candidate hits.[m
[32m+[m[32m5. Optionally re-rank them using Smith-Waterman score.[m
[32m+[m
[32m+[m[32mThis makes the project more than just a pairwise alignment implementation; it becomes a basic search pipeline.[m
[32m+[m
[32m+[m[32m## Benchmarks[m
[32m+[m
[32m+[m[32mBenchmarking results are documented in:[m
[32m+[m
[32m+[m[32m```text[m
[32m+[m[32mbenchmarks/BENCHMARKS.md[m
[32m+[m[32m```[m
[32m+[m
[32m+[m[32mCurrent benchmarks evaluate Smith-Waterman runtime scaling on ASTRAL/SCOPe protein sequence datasets.[m
[32m+[m
[32m+[m[32mThe benchmark report currently compares:[m
[32m+[m
[32m+[m[32m- score-only Smith-Waterman computation[m
[32m+[m[32m- exhaustive alignment reconstruction[m
[32m+[m[32m- runtime scaling across increasing dataset sizes[m
[32m+[m[32m- total residue count as an explanation for runtime growth[m
[32m+[m
[32m+[m[32mThe benchmarked datasets include chunks of:[m
[32m+[m
[32m+[m[32m- 10 sequences[m
[32m+[m[32m- 100 sequences[m
[32m+[m[32m- 1000 sequences[m
[32m+[m[32m- 10000 sequences[m
[32m+[m
[32m+[m[32mThe benchmark report also documents current limitations, including:[m
[32m+[m
[32m+[m[32m- linear gap model only[m
[32m+[m[32m- simple scoring system[m
[32m+[m[32m- no affine gap penalties yet[m
[32m+[m[32m- exact Smith-Waterman benchmark only[m
[32m+[m[32m- no current k-mer search benchmark yet[m
[32m+[m
[32m+[m[32m## Running Tests[m
[32m+[m
[32m+[m[32mRun all tests:[m
[32m+[m
[32m+[m[32m```bash[m
[32m+[m[32mpytest[m
[32m+[m[32m```[m
[32m+[m
[32m+[m[32mRun only search pipeline tests:[m
[32m+[m
[32m+[m[32m```bash[m
[32m+[m[32mpytest tests/pipelines/[m
[32m+[m[32m```[m
[32m+[m
[32m+[m[32mRun only benchmark tests:[m
[32m+[m
[32m+[m[32m```bash[m
[32m+[m[32mpytest tests/benchmarks/[m
[32m+[m[32m```[m
[32m+[m
[32m+[m[32mThe tests currently cover:[m
[32m+[m
[32m+[m[32m- DNA utility functions[m
[32m+[m[32m- database normalization[m
[32m+[m[32m- Needleman-Wunsch alignment[m
[32m+[m[32m- Smith-Waterman alignment[m
[32m+[m[32m- k-mer search[m
[32m+[m[32m- search refinement[m
[32m+[m[32m- full search pipeline behavior[m
[32m+[m[32m- benchmark smoke tests[m
[32m+[m[32m- benchmark residue-count checks[m
[32m+[m
[32m+[m[32m## Examples[m
[32m+[m
[32m+[m[32mA simple search example is provided in:[m
[32m+[m
[32m+[m[32m```text[m
[32m+[m[32mexamples/search_demo.py[m
[32m+[m[32m```[m
[32m+[m
[32m+[m[32mRun it with:[m
[32m+[m
[32m+[m[32m```bash[m
[32m+[m[32mpython examples/search_demo.py[m
[32m+[m[32m```[m
[32m+[m
[32m+[m[32m## Current Strengths[m
[32m+[m
[32m+[m[32mThe strongest parts of the project are:[m
[32m+[m
[32m+[m[32m- Clear movement from isolated algorithms toward a search pipeline[m
[32m+[m[32m- Basic package organization[m
[32m+[m[32m- Tests for search, refinement, database normalization, and benchmark behavior[m
[32m+[m[32m- Benchmark documentation instead of only toy examples[m
[32m+[m[32m- Honest educational scope[m
[32m+[m[32m- Beginning of reproducibility through dataset chunks and benchmark reports[m
[32m+[m
[32m+[m[32mThe project is especially useful for learning how sequence database search can be built from smaller algorithmic pieces.[m
[32m+[m
[32m+[m[32m## Current Limitations[m
[32m+[m
[32m+[m[32mThis project is still early-stage. Important limitations include:[m
[32m+[m
[32m+[m[32m- No BLOSUM or PAM substitution matrix support yet[m
[32m+[m[32m- No affine gap penalties yet[m
[32m+[m[32m- No statistical significance estimates such as E-values or bit scores[m
[32m+[m[32m- No indexed k-mer search yet[m
[32m+[m[32m- Current k-mer search scans database sequences directly[m
[32m+[m[32m- Alignment result objects are still minimal[m
[32m+[m[32m- No stable command-line interface yet[m
[32m+[m[32m- No biological case study has been completed yet[m
[32m+[m[32m- Not intended for production biological analysis[m
[32m+[m
[32m+[m[32mThese limitations are intentional development targets, not hidden assumptions.[m
[32m+[m
[32m+[m[32m## Roadmap[m
[32m+[m
[32m+[m[32mPlanned development stages:[m
[32m+[m
[32m+[m[32m### 1. Improve Documentation and Reproducibility[m
[32m+[m
[32m+[m[32m- Expand README documentation[m
[32m+[m[32m- Keep benchmark results synchronized with benchmark scripts[m
[32m+[m[32m- Add clearer instructions for running examples and benchmarks[m
[32m+[m
[32m+[m[32m### 2. Improve Alignment Output[m
[32m+[m
[32m+[m[32m- Return structured alignment result objects[m
[32m+[m[32m- Include alignment score[m
[32m+[m[32m- Include aligned sequences[m
[32m+[m[32m- Include identity percentage[m
[32m+[m[32m- Include gap counts[m
[32m+[m[32m- Include local alignment start/end positions where possible[m
[32m+[m
[32m+[m[32m### 3. Add Biologically Realistic Scoring[m
[32m+[m
[32m+[m[32m- Add BLOSUM62 support for protein alignment[m
[32m+[m[32m- Add support for additional scoring matrices later[m
[32m+[m[32m- Add affine gap penalties with separate gap opening and gap extension costs[m
[32m+[m
[32m+[m[32m### 4. Benchmark Search Strategies[m
[32m+[m
[32m+[m[32mCompare:[m
[32m+[m
[32m+[m[32m- exhaustive Smith-Waterman search[m
[32m+[m[32m- k-mer search only[m
[32m+[m[32m- k-mer search with Smith-Waterman refinement[m
[32m+[m
[32m+[m[32mMeasure:[m
[32m+[m
[32m+[m[32m- runtime[m
[32m+[m[32m- dataset size[m
[32m+[m[32m- total residue count[m
[32m+[m[32m- top-hit recovery[m
[32m+[m[32m- speed/sensitivity tradeoffs[m
[32m+[m
[32m+[m[32m### 5. Add Indexed Search[m
[32m+[m
[32m+[m[32mBuild an inverted k-mer index so that database k-mers do not need to be regenerated for every query.[m
[32m+[m
[32m+[m[32mExample idea:[m
[32m+[m
[32m+[m[32m```python[m
[32m+[m[32m{[m
[32m+[m[32m    "ATG": ["seq1", "seq5", "seq9"],[m
[32m+[m[32m    "TGC": ["seq1", "seq2"][m
[32m+[m[32m}[m
[32m+[m[32m```[m
[32m+[m
[32m+[m[32mThis would make repeated searches more scalable.[m
[32m+[m
[32m+[m[32m### 6. Add a Biological Case Study[m
[32m+[m
[32m+[m[32mUse the toolkit on a real protein-family dataset.[m
[32m+[m
[32m+[m[32mPossible case study direction:[m
[32m+[m
[32m+[m[32m- collect related protein sequences[m
[32m+[m[32m- search against a local FASTA database[m
[32m+[m[32m- refine hits with Smith-Waterman[m
[32m+[m[32m- compare conserved regions[m
[32m+[m[32m- discuss biological interpretation and limitations[m
[32m+[m
[32m+[m[32m## Educational Scope[m
[32m+[m
[32m+[m[32mThis repository is meant to show the internal logic behind common bioinformatics sequence-analysis tasks.[m
[32m+[m
[32m+[m[32mIt is not currently designed for clinical, diagnostic, or production research use.[m
[32m+[m
[32m+[m[32mFor serious biological analysis, established tools such as BLAST, HMMER, EMBOSS, MAFFT, MUSCLE, Clustal Omega, and Biopython should be used. This project is mainly for learning, experimentation, and building algorithmic understanding.[m
[32m+[m
[32m+[m[32m## Dependencies[m
[32m+[m
[32m+[m[32mCurrent dependencies are listed in:[m
[32m+[m
[32m+[m[32m```text[m
[32m+[m[32mrequirements.txt[m
[32m+[m[32m```[m
[32m+[m
[32m+[m[32mMain dependencies:[m
[32m+[m
[32m+[m[32m- Biopython[m
[32m+[m[32m- pytest[m
[32m+[m[32m- pandas[m
[32m+[m
[32m+[m[32m## License[m
[32m+[m
[32m+[m[32mNo license has been specified yet.[m
\ No newline at end of file[m
