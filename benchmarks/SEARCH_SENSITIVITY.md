# Search Sensitivity Benchmark Report

## Purpose

This report evaluates the search-sensitivity behavior of the current `bioseq-toolkit` search pipeline.

The comparison is between:

1. exhaustive Smith-Waterman search,
2. k-mer-only search,
3. k-mer search followed by Smith-Waterman refinement.

The goal is to test whether the heuristic k-mer stage can preserve the high-confidence exhaustive Smith-Waterman hits while reducing the number of database sequences passed into expensive local alignment.

This is an internal algorithmic benchmark. Exhaustive Smith-Waterman is used as the reference ranking, not as biological ground truth.

---

## Benchmark Setup

### Dataset family

The benchmarks use ASTRAL/SCOPe protein sequence subsets from:

```text
data/benchmark_sequences/
```

The main tested subsets are:

```text
astral_100.fasta
astral_1000.fasta
astral_10000.fasta
```

### Query

The same protein query was used across benchmark runs:

```text
ttltnpqkaairsswskfmdngvsngqgfymdlfkahpetltpfkslfggltlaqlqdnpkmkaqslvfcngmssfvdhlddndmlvvliqkmaklhnnrgirasdlrtaydilihymedhnhmvggakdawevfvgficktlgdymkels
```

### Scoring

```text
Matrix: BLOSUM62
Gap penalty: -4
Gap model: linear
```

### Reference grouping

Exhaustive Smith-Waterman scores were used to build reference tiers automatically.

For this query, the high-confidence set is:

```text
Tier 1 exact/near-exact hits: 15 sequences
Tier 2 high-scoring hits:      5 sequences
Tier 1 + Tier 2 total:        20 high-confidence sequences
```

The main sensitivity metric is recovery of these 20 high-confidence hits.

---

## Main Result

Across the benchmarked datasets, `k = 3` with threshold `1` gave the best current tradeoff.

|     Dataset   | Database Size | k | Threshold | Candidates | Candidate Fraction | SW Recall@20 | Tier 1+2 Recall | Background Returned |
|:-------------:|:-------------:|:-:|:---------:|:----------:|:------------------:|:------------:|:---------------:|:-------------------:|
| `astral_100`  |      100      | 3 |     1     |  20 / 100  |        20.0%       |    20 / 20   |     20 / 20     |           0         |
| `astral_1000` |      1000     | 3 |     1     |  20 / 1000 |        2.0%        |    20 / 20   |     20 / 20     |           0         |
| `astral_10000`|      10000    | 3 |     1     |  20 / 10000|        0.2%        |    20 / 20   |     20 / 20     |           0         |

The strongest observation is that the candidate count stayed at 20 while the database grew from 100 to 10000 sequences.

---

## `astral_100` Result

### Best setting

```text
Dataset: data/benchmark_sequences/astral_100.fasta
Database size: 100 sequences
k: 3
Threshold: 1
Matrix: BLOSUM62
Gap penalty: -4
```

### Recovery

|         Method        | Candidates | Recall@5 | Recall@10 | Recall@20 | Tier 1+2 Recall | Background Returned |
|:---------------------:|:----------:|:--------:|:---------:|:---------:|:---------------:|:-------------------:|
|       K-mer only      |  20 / 100  |   5 / 5  |  10 / 10  |  20 / 20  |     20 / 20     |           0         |
| K-mer + SW refinement |  20 / 100  |   5 / 5  |  10 / 10  |  20 / 20  |     20 / 20     |           0         |

### Runtime

### Runtime

### Runtime

|         Method        |   Runtime  | Speedup vs exhaustive SW |
|:---------------------:|:----------:|:------------------------:|
| Exhaustive SW         | 5.029279 s |           1.00x          |
| K-mer only            | 0.003196 s |           1573.47x       |
| K-mer + SW refinement | 1.060561 s |           4.74x          |

### Interpretation

For `astral_100`, k-mer-only search recovered the same top 20 high-confidence hits as exhaustive SW while reducing the database from 100 sequences to 20 candidates.

---

## Parameter Sweep — `astral_100`

A sweep was run across:

```text
K_VALUES = [2, 3, 4, 5]
THRESHOLDS = [1, 2, 3, 5, 10]
```

### Summary by k-value

| k |                 Candidate Behavior                  |         Recall Behavior         |             Interpretation             |
|:-:|:---------------------------------------------------:|:-------------------------------:|:--------------------------------------:|
| 2 | Returned 86 / 100 candidates and 39 background hits | 20 / 20 high-confidence recall  | Sensitive but too broad.               |
| 3 | Returned 20 / 100 candidates and 0 background hits  | 20 / 20 high-confidence recall  | Best tradeoff.                         |
| 4 | Returned 15 / 100 candidates                        | 15 / 20 high-confidence recall  | Too strict; loses Tier 2 hits.         |
| 5 | Returned 15 / 100 candidates                        | 15 / 20 high-confidence recall  | Too strict; same issue as k = 4.       |

### Sweep interpretation

`k = 2` preserved recall but returned too much of the database. `k = 4` and `k = 5` were stricter but lost the five Tier 2 high-scoring hits. `k = 3` preserved all high-confidence hits while returning no background hits.

Threshold had little effect in this benchmark compared with k-mer size.

---

## Larger Dataset Check — `astral_1000`

The `astral_1000` run was used as an intermediate stress test between the 100-sequence and 10000-sequence datasets.

The important verified pattern was:

```text
k = 3
threshold = 1
candidates = 20 / 1000
candidate fraction = 2.0%
SW Recall@20 = 20 / 20
Tier 1+2 recall = 20 / 20
background returned = 0
```

The main qualitative result was the same as `astral_100`: `k = 3` preserved the high-confidence SW hits while reducing the candidate set strongly.

A separate observation from the larger sweep was that `k = 2` became much slower during refinement because it retained too many candidates. This made k-mer + SW refinement approach exhaustive-search behavior for permissive k-mer settings.

---

## Large Dataset Result — `astral_10000`

After `k = 3` was selected from the smaller runs, the benchmark was run on `astral_10000` using the selected setting only.

### Settings

```text
Dataset: data/benchmark_sequences/astral_10000.fasta
Database size: 10000 sequences
k: 3
Threshold: 1
Matrix: BLOSUM62
Gap penalty: -4
```

### Recovery

|         Method          |  Candidates  | Recall@5 | Recall@10 | Recall@20 | Tier 1+2 Recall | Background Returned |
|:-----------------------:|:------------:|:--------:|:---------:|:---------:|:---------------:|:-------------------:|
|       K-mer only        | 20 / 10000   |  5 / 5   |  10 / 10  |  20 / 20  |     20 / 20     |          0          |
| K-mer + SW refinement   | 20 / 10000   |  5 / 5   |  10 / 10  |  20 / 20  |     20 / 20     |          0          |

### Runtime

|         Method          |    Runtime    | Speedup vs exhaustive SW |
|:-----------------------:|:-------------:|:------------------------:|
|      Exhaustive SW      | 639.680913 s  |          1.00x           |
|       K-mer only        |  0.527722 s   |        1212.16x          |
| K-mer + SW refinement   |  1.724182 s   |         371.01x          |

### Interpretation

The `astral_10000` run is the strongest result so far. With `k = 3` and threshold `1`, the k-mer stage returned only 20 candidates out of 10000 sequences while preserving full recovery of the exhaustive SW top 20 and the full Tier 1+2 high-confidence set.

The exhaustive SW runtime was reused from the cached exhaustive-SW reference-generation run. The cache stores per-sequence SW scores and lengths, while the benchmark summary CSV stores recall, candidate count, runtime, and speedup.

---

## Cross-Dataset Summary

Selected setting:

```text
k = 3
threshold = 1
matrix = BLOSUM62
gap penalty = -4
```

|     Dataset      | Database Size |  Candidates  | Candidate Fraction | SW Recall@20 | Tier 1+2 Recall | Background Returned |
|:----------------:|:-------------:|:------------:|:------------------:|:------------:|:---------------:|:-------------------:|
|  `astral_100`    |      100      |   20 / 100   |       20.0%        |   20 / 20    |     20 / 20     |          0          |
|  `astral_1000`   |     1000      |  20 / 1000   |        2.0%        |   20 / 20    |     20 / 20     |          0          |
| `astral_10000`   |     10000     | 20 / 10000   |        0.2%        |   20 / 20    |     20 / 20     |          0          |

|     Dataset      | Exhaustive SW Runtime |   K-mer Runtime  | Refined Total Runtime |
|:----------------:|:---------------------:|:----------------:|:---------------------:|
|  `astral_100`    |      5.029279 s       |   0.003196 s     |      1.060561 s       |
|  `astral_1000`   |   not recorded here   | not recorded here|   not recorded here  |
| `astral_10000`   |     639.680913 s      |   0.527722 s     |      1.724182 s       |

---

## Main Conclusion

For the current query and ASTRAL benchmark subsets, `k = 3` with threshold `1` is the best current k-mer setting.

It preserved:

```text
SW Recall@5:     5 / 5
SW Recall@10:   10 / 10
SW Recall@20:   20 / 20
Tier 1+2 recall: 20 / 20
Background returned: 0
```

while sharply reducing candidate counts.

The main parameter pattern is:

```text
k = 2  -> too permissive; high recall but many candidates/background hits.
k = 3  -> best tradeoff; full high-confidence recall with low candidate count.
k = 4  -> too strict; loses Tier 2 high-scoring hits.
k = 5  -> too strict; same issue as k = 4.
```

---

## Refinement Interpretation

SW refinement is useful in principle because it can re-rank candidates using a stronger local-alignment score. However, in the current benchmark, k-mer-only search already recovered the same top 20 high-confidence hits as exhaustive SW for the selected `k = 3` setting.

Therefore, in these runs:

```text
K-mer-only search recovered the high-confidence set.
K-mer + SW refinement remained much faster than exhaustive SW, but did not improve recall in this specific benchmark.
```

The value of refinement should be tested later on harder queries or more diverse datasets where k-mer-only ranking begins to fail.

---

## Seed Extension Decision

Seed extension is not urgent yet.

The current simple k-mer method already achieved full high-confidence recovery across the tested datasets with `k = 3`.

Seed extension becomes meaningful only after identifying a failure case, such as:

```text
k-mer-only misses Tier 2 hits,
k-mer-only returns too many background hits,
k-mer-only ranking differs from exhaustive SW ranking,
k-mer-only performs badly on a harder dataset or query.
```

The better next step is to preserve this benchmark as the baseline, then test harder biological cases or implement indexed search.

---

## Files Produced

Recommended report path:

```text
benchmarks/SEARCH_SENSITIVITY.md
```

Recommended raw result files:

```text
benchmarks/search_sensitivity_results_astral_10000_k3_threshold1.csv
benchmarks/cached_exhaustive_sw_results/astral_10000_BLOSUM62_gap-4.csv
```

---

## Limitations

This benchmark should not be overgeneralized.

Important limitations:

1. The benchmark uses one query sequence.
2. Exhaustive SW is used as an internal reference ranking, not biological ground truth.
3. The high-confidence tiers are derived from SW score thresholds, not curated family labels.
4. The current k-mer method is still scan-based rather than indexed.
5. The current implementation uses a linear gap penalty.
6. The strong result may partly reflect redundancy or close homologs in the ASTRAL subset.

---

## Next Steps

Recommended immediate next steps:

1. Commit the sensitivity benchmark script and final report.
2. Keep the cached exhaustive SW CSV for `astral_10000`.
3. Do not implement seed extension yet.
4. Next improvement should likely be indexed k-mer search or a real biological family case study.
