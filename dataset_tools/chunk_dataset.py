from pathlib import Path

from bioseq.fasta_io import read_fasta_records


PROJECT_ROOT = Path(__file__).resolve().parent.parent #Always goes to root, no matter w

OUTPUT_DIR = (
    PROJECT_ROOT
    / "data"
    / "benchmark_sequences"
)


def dataset_to_chunks(dataset, *chunk_lengths):
    """
    Split a FASTA dataset into multiple benchmark
    subsets and write them as FASTA files.

    Parameters
    ----------
    dataset: str
        Dataset path
    *chunk_lengths : int
        One or more chunk sizes specifying how many
        sequences should be written into each output dataset.

    
    """
    dataset_sequences = read_fasta_records(dataset)

    dataset_name = Path(dataset).stem

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    for chunk_size in chunk_lengths:

        temp_chunk = dataset_sequences[:chunk_size]

        output_file = (
            OUTPUT_DIR
            / f"{dataset_name}_{chunk_size}.fasta"
        )

        with open(output_file, "w") as file:

            for record in temp_chunk:

                file.write(
                    record["header"] + "\n"
                )

                file.write(
                    record["sequence"] + "\n"
                )