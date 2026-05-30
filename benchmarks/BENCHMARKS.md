# Benchmarking Results

## Objective

Evaluate how Smith-Waterman search scales with increasing database size and compare the computational cost of exhaustive alignment reconstruction versus score-only computation.

## Experimental Setup

Dataset:

* ASTRAL/SCOPe protein sequences

Query:

* Single protein sequence (~142 residues)

Algorithm:

* Smith-Waterman local alignment

Scoring:

* Linear gap penalties
* Simple scoring system

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

## Observations

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
| 1000    |        142,343 |
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
* Current benchmark measures exact search only
* Hardware specifications not controlled/reported

## Future Benchmarks

* Single traceback Smith-Waterman
* k-mer filtering
* k-mer + refinement pipeline
* Indexed search
* Affine gap penalties
* Runtime per residue analysis
* Sensitivity/speed tradeoff experiments
