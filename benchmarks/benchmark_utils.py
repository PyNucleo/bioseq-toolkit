from bioseq.fasta_io import read_fasta_records

def get_average_runtime(benchmark_method, iterations, dataset):
    """
    Run a benchmark function multiple times and return
    the average runtime.

    Parameters
    ----------
    benchmark_method : callable
        Benchmark function to execute.
    iterations : int
        Number of repeated runs.
    dataset : str
        Dataset path passed into benchmark_method.

    Returns
    -------
    float
        Mean runtime across all benchmark iterations.
    """

    total_runtime = 0

    for _ in range(iterations):
        runtime = benchmark_method(dataset)
        total_runtime += runtime

    return total_runtime / iterations

def get_total_dataset_residues(dataset):
    dataset_records = read_fasta_records(dataset)

    temp_s = 0

    for seq in dataset_records:
        temp_s += len(seq["sequence"])
    
    return temp_s
