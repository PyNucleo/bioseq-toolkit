from benchmarks.benchmark_alignment import benchmark_smith_waterman_scores_run_time, benchmark_smith_waterman_alignment_run_time
from benchmarks.benchmark_utils import get_total_dataset_residues

def test_score_only_benchmark(tmp_path):

    fasta = tmp_path / "test.fasta" 

    fasta.write_text(
        ">sq1\n"
        "ATGCHA\n"
        ">sq2\n"
        "YFWGGG\n"
        ">sq3\n"
        "GCHIKL\n"
        ">sq4\n"
        "ATGCHAFTYUIAC"
    )

    runtime = benchmark_smith_waterman_scores_run_time(fasta)

    assert isinstance(runtime, float)
    assert runtime >= 0

def test_alignment_reconstruction_benchmark(tmp_path):

    fasta = tmp_path / "tests.fasta"

    fasta.write_text(
        ">sq1\n"
        "ATGCHA\n"
        ">sq2\n"
        "YFWGGG\n"
        ">sq3\n"
        "GCHIKL\n"
        ">sq4\n"
        "ATGCHAFTYUIAC"
    )

    runtime = benchmark_smith_waterman_alignment_run_time(fasta)

    assert isinstance(runtime, float)
    assert runtime >= 0

def test_total_dataset_residues():

    dataset_10 = "data/benchmark_sequences/astral_10.fasta"
    dataset_100 = "data/benchmark_sequences/astral_100.fasta"
    dataset_1000 = "data/benchmark_sequences/astral_1000.fasta"
    dataset_10000 = "data/benchmark_sequences/astral_10000.fasta"

    residues_in_dataset_10 = get_total_dataset_residues(dataset_10)
    residues_in_dataset_100 = get_total_dataset_residues(dataset_100)
    residues_in_dataset_1000 = get_total_dataset_residues(dataset_1000)
    residues_in_dataset_10000 = get_total_dataset_residues(dataset_10000)

    assert residues_in_dataset_10 == 1401
    assert residues_in_dataset_100 == 14274
    assert residues_in_dataset_1000 == 142363
    assert residues_in_dataset_10000 == 2471135




