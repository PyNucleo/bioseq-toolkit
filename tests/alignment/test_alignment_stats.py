import pytest

from bioseq.alignment.alignment_stats import (
    get_alignment_stats,
    get_identity,
    get_matches_mismatches_gaps,
)


def test_alignment_stats_for_gapped_alignment():
    aligned_1 = "ATGCG-"
    aligned_2 = "AT-CGA"

    stats = get_alignment_stats(aligned_1, aligned_2)

    assert stats == {
        "alignment_length": 6,
        "matches": 4,
        "mismatches": 0,
        "gaps": 2,
        "gap_columns": 2,
        "identity": 0.667,
        "identity_excluding_gaps": 1.0,
        "similarity": None,
    }


def test_alignment_stats_for_mismatches_without_gaps():
    aligned_1 = "ATGCG"
    aligned_2 = "ATCGA"

    stats = get_alignment_stats(aligned_1, aligned_2)

    assert stats["alignment_length"] == 5
    assert stats["matches"] == 2
    assert stats["mismatches"] == 3
    assert stats["gaps"] == 0
    assert stats["gap_columns"] == 0
    assert stats["identity"] == 0.4
    assert stats["identity_excluding_gaps"] == 0.4
    assert stats["similarity"] is None


def test_alignment_stats_for_perfect_alignment():
    aligned_1 = "GATTACA"
    aligned_2 = "GATTACA"

    stats = get_alignment_stats(aligned_1, aligned_2)

    assert stats["alignment_length"] == 7
    assert stats["matches"] == 7
    assert stats["mismatches"] == 0
    assert stats["gaps"] == 0
    assert stats["gap_columns"] == 0
    assert stats["identity"] == 1.0
    assert stats["identity_excluding_gaps"] == 1.0


def test_get_identity_including_and_excluding_gaps():
    aligned_1 = "ATGCG-"
    aligned_2 = "AT-CGA"

    assert get_identity(aligned_1, aligned_2, including_gaps=True) == 0.667
    assert get_identity(aligned_1, aligned_2, including_gaps=False) == 1.0


def test_get_matches_mismatches_gaps():
    aligned_1 = "ATGCG-"
    aligned_2 = "AT-CGA"

    matches, mismatches, gaps, gap_columns = get_matches_mismatches_gaps(
        aligned_1,
        aligned_2,
    )

    assert matches == 4
    assert mismatches == 0
    assert gaps == 2
    assert gap_columns == 2


def test_alignment_stats_rejects_unequal_lengths():
    with pytest.raises(ValueError):
        get_alignment_stats("ATGC", "ATGCG")