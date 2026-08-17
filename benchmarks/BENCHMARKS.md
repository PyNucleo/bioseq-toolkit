# Benchmarking Results

## Objective

Evaluate how Smith-Waterman search scales with increasing database size and compare the computational cost of exhaustive alignment reconstruction, score-only computation, k-mer search, and k-mer search with Smith-Waterman refinement.

## Recorded Setup

Dataset:

* ASTRAL/SCOPe protein sequences

Query:

* Single hard-coded protein sequence (151 residues)

Algorithm:

* Smith-Waterman local alignment

Scoring:

* Linear gap penalties
* Alignment reconstruction and score-only timing: match `2`, mismatch `-1`, gap `-2`
* Refined public-search timing: `search()` defaults match `1`, mismatch `-1`, gap `-2`
* The separate sensitivity report uses BLOSUM62 with gap `-4`

Hardware:
- Machine: Lenovo IdeaPad Slim 3i
- CPU: 13th Gen Intel Core i7-13620H
- Base Clock: 2.40 GHz
- Cores / Threads: 10 cores / 16 logical processors
- RAM: 16.0 GB DDR5 (15.7 GB usable)
- Memory Speed: 4800 MT/s

Benchmarks were performed on a consumer-grade laptop system under normal operating conditions; background processes were not strictly controlled.

Software Environment:
- Operating System: Microsoft Windows 11 Version 24H2 (Build 26100.2894)
- Python Version: 3.14.5

Runs:

* chunk_10: 100 iterations
* chunk_100: 50 iterations
* chunk_1000: 25 iterations
* chunk_10000: 10 iterations

Benchmark Modes:

* Exhaustive Smith-Waterman alignment reconstruction
* Score-only Smith-Waterman computation
* K-mer-only search
* K-mer search with Smith-Waterman refinement

The numbers below are existing recorded benchmark results from the scripts and environment described above. They were not regenerated as part of this documentation update.

`benchmark_alignment.py` and `benchmark_search.py` contain executable timing
functions; the `run_*_benchmarks.py` modules execute them. `benchmarks/main.py`
only prints static historical values and does not rerun a benchmark.

## Dataset Statistics

| Dataset | Total Residues | Residue Increase |
| ------- | -------------: | ---------------: |
| 10      |          1,401 |                — |
| 100     |         14,274 |           10.19× |
| 1000    |        142,363 |            9.97× |
| 10000   |      2,471,135 |           17.36× |

## Runtime Results — Score Only

| Dataset | Runtime (s) | Runtime Increase |
| ------- | ----------: | ---------------: |
| 10      |       0.078 |                — |
| 100     |       0.930 |            11.9× |
| 1000    |      10.765 |            11.6× |
| 10000   |     137.183 |            12.7× |

## Runtime Results — Exhaustive Alignment Reconstruction

| Dataset | Runtime (s) | Runtime Increase |
| ------- | ----------: | ---------------: |
| 10      |       0.115 |                — |
| 100     |       1.151 |            10.0× |
| 1000    |      11.825 |            10.3× |
| 10000   |     218.914 |            18.5× |

## Alignment vs Score-Only Comparison

| Dataset | Alignment (s) | Score-Only (s) | Slowdown |
| ------- | ------------: | -------------: | -------: |
| 10      |         0.115 |          0.078 |    1.47× |
| 100     |         1.151 |          0.930 |    1.24× |
| 1000    |        11.825 |         10.765 |    1.10× |
| 10000   |       218.914 |        137.183 |    1.60× |

## Alignment Runtime Observations

* Runtime scaled approximately linearly from 10 → 1000 sequences for both benchmark modes.
* Runtime increase between 1000 and 10000 sequences was larger than expected from sequence count alone.
* Initially, exhaustive traceback was suspected to be the dominant contributor because enumerating many optimal alignments can significantly increase computational cost.
* However, the same scaling anomaly persisted in score-only computation, indicating traceback was not the primary explanation.

## Investigation of Runtime Jump

The primary explanation was total residue count rather than sequence count.

| Dataset | Total Residues |
| ------- | -------------: |
| 10      |          1,401 |
| 100     |         14,274 |
| 1000    |        142,363 |
| 10000   |      2,471,135 |

Residue count increased by approximately 17.4× between the 1000 and 10000 datasets rather than the expected ~10× increase from sequence count alone.

Since Smith-Waterman complexity is dominated by:

O(query_length × target_length)

runtime tracked residue growth more closely than sequence count.

Residue count for each dataset was computed by iterating through all sequences and summing sequence lengths.

## Runtime Scaling Visualization

![Runtime scaling](figures/sw_score_runtime_vs_residues.png)

The visualization suggests runtime scales approximately proportionally with residue count, supporting the interpretation that total search space rather than sequence count dominated computational cost.

## Limitations

* Linear gap model only
* Protein lengths vary substantially between datasets
* No affine gap penalties
* Current benchmarks emphasize runtime, not biological sensitivity
* Hardware specifications are reported, but hardware and background processes were not strictly controlled

## Future Benchmarks

* Regular-versus-indexed multi-search performance comparison
* Affine gap penalties
* Runtime per residue analysis
* Sensitivity/speed tradeoff experiments
* Seed-extension search variants

## Runtime Results - K-mer Search Only Using Varying k Sizes and Thresholds 

These historical measurements call the scan-based `kmer_search()` function.
Current source also supports indexed multi-search, but this report contains no
regular-versus-indexed measurement and makes no indexed-speed claim.

| Dataset      |  k | Threshold | Runtime (s) |
| ------------ | -: | --------: | ----------: |
| astral_10    |  3 |         1 |    0.000608 |
| astral_10    |  3 |         3 |    0.000538 |
| astral_10    |  3 |         5 |    0.000546 |
| astral_10    |  4 |         2 |    0.000511 |
| astral_10    |  4 |         3 |    0.000511 |
| astral_100   |  3 |         1 |    0.003706 |
| astral_100   |  3 |         3 |    0.003846 |
| astral_100   |  3 |         5 |    0.003831 |
| astral_100   |  4 |         2 |    0.003782 |
| astral_100   |  4 |         3 |    0.004647 |
| astral_1000  |  3 |         1 |    0.033992 |
| astral_1000  |  3 |         3 |    0.034162 |
| astral_1000  |  3 |         5 |    0.033632 |
| astral_1000  |  4 |         2 |    0.031665 |
| astral_1000  |  4 |         3 |    0.031457 |
| astral_10000 |  3 |         1 |    0.520231 |
| astral_10000 |  3 |         3 |    0.520493 |
| astral_10000 |  3 |         5 |    0.522213 |
| astral_10000 |  4 |         2 |    0.495486 |
| astral_10000 |  4 |         3 |    0.495018 |

## K-mer Search Observations

K-mer-only search was substantially faster than exact Smith-Waterman score-only search. On astral_10000, exact SW score-only took about 137.18 s, while k-mer-only search took about 0.50–0.52 s depending on k and threshold. 

Changing the threshold had little effect on runtime because the current implementation still scans all database sequences before filtering. Increasing k from 3 to 4 slightly reduced runtime on the largest dataset, likely because fewer k-mers are generated per sequence, which further supports the impact of kmer-filtering on larger datasets when it comes to cutting down runtime.

## Runtime Results — K-mer + Smith-Waterman Refinement

Parameters:
- k = 3
- threshold = 3
- top_n_hits = 10
- simple refinement scoring = match 1, mismatch -1, linear gap -2

| Dataset | Exact SW Score-Only (s) | K-mer + SW Refinement (s) |  Speedup  |
|---------|-------------------------|---------------------------|-----------|
|    10   |          0.078          |           0.036           |   2.18×   |
|    100  |          0.930          |           0.092           |   10.13×  |
|    1000 |          10.765         |           0.126           |   85.21×  |
|    10000|          137.183        |           0.673           |   203.73× |

## Interpretation

K-mer + Smith-Waterman refinement was much faster than exhaustive Smith-Waterman scoring because Smith-Waterman was only applied to the top candidate hits rather than every database sequence.

This does not prove equal sensitivity or biological accuracy. It only demonstrates runtime improvement for the current query, datasets, and parameters.

## Search Sensitivity Benchmark

A separate sensitivity benchmark was added in [`SEARCH_SENSITIVITY.md`](SEARCH_SENSITIVITY.md).

That report compares exhaustive Smith-Waterman, k-mer-only search, and k-mer search with Smith-Waterman refinement by measuring top-hit recovery, high-confidence hit recovery, candidate-set reduction, background returned, and runtime.

The sensitivity benchmark is kept separate from this file because this document primarily records runtime-scaling results, while `SEARCH_SENSITIVITY.md` focuses on recovery and speed-sensitivity tradeoffs.
