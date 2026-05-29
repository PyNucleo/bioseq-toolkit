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