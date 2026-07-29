import pytest

from database.database_utils import normalize_database
from database.sequence_database import SequenceDatabase


def test_normalize_database():
    raw_sequences = [
        "ATGCGT",
        "ATGCGT",
        "GGGGGG"
    ]

    normalized = normalize_database(raw_sequences)
    records = normalized.get_sequences()

    assert isinstance(normalized, SequenceDatabase)

    assert records == [
        {"id": "id1", "sequence": "ATGCGT"},
        {"id": "id2", "sequence": "ATGCGT"},
        {"id": "id3", "sequence": "GGGGGG"},
    ]

    assert records[0]["sequence"] == records[1]["sequence"]
    assert records[0]["id"] != records[1]["id"]

    existing_db = SequenceDatabase([
        {"id": "custom", "sequence": "AAAA"}
    ])

    assert normalize_database(existing_db) is existing_db

    with pytest.raises(TypeError):
        normalize_database(123)
    with pytest.raises(TypeError):
        normalize_database({"seq1": "ATGC"})


def test_normalize_database_rejects_duplicate_ids_with_all_positions():
    database = SequenceDatabase([
        {"id": "alpha", "sequence": "AAAA"},
        {"id": "beta", "sequence": "CCCC"},
        {"id": "alpha", "sequence": "GGGG"},
        {"id": "beta", "sequence": "TTTT"},
        {"id": "beta", "sequence": "ACGT"},
    ])

    with pytest.raises(ValueError) as error:
        normalize_database(database)

    message = str(error.value)
    assert "Duplicate database IDs detected:" in message
    assert "'alpha' appears in records [1, 3]" in message
    assert "'beta' appears in records [2, 4, 5]" in message


def test_normalize_database_rejects_duplicate_ids_from_fasta(tmp_path):
    fasta = tmp_path / "duplicate_ids.fasta"
    fasta.write_text(
        ">duplicate first record\n"
        "AAAA\n"
        ">unique\n"
        "CCCC\n"
        ">duplicate second record\n"
        "GGGG\n"
    )

    with pytest.raises(ValueError, match="'duplicate' appears in records \\[1, 3\\]"):
        normalize_database(str(fasta))


def test_normalize_database_accepts_unique_ids_with_repeated_sequences():
    database = SequenceDatabase([
        {"id": "first", "sequence": "AAAA"},
        {"id": "second", "sequence": "AAAA"},
    ])

    assert normalize_database(database) is database
