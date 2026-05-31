import pytest

from bioseq.alignment.smith_waterman import local_alignment


def assert_basic_structured_result(result, s1, s2, score, best_positions):
    assert result["algorithm"] == "Smith-Waterman"
    assert result["mode"] == "local"
    assert result["sequence_1"] == s1
    assert result["sequence_2"] == s2
    assert result["score"] == score
    assert result["best_positions"] == best_positions

    assert result["scoring"]["matrix"] is None
    assert result["scoring"]["gap_model"] == "linear"

    assert isinstance(result["num_alignments"], int)
    assert isinstance(result["alignments"], list)
    assert result["num_alignments"] == len(result["alignments"])


def test_local_alignment_identical_sequences_structured():
    result = local_alignment(
        "ATGC",
        "ATGC",
        match=2,
        mismatch=-1,
        gap_penalty=-2,
        return_all=False,
        structured=True,
    )

    assert_basic_structured_result(
        result,
        s1="ATGC",
        s2="ATGC",
        score=8,
        best_positions=[(4, 4)],
    )

    assert result["scoring"] == {
        "match": 2,
        "mismatch": -1,
        "gap": -2,
        "matrix": None,
        "gap_model": "linear",
    }

    assert result["num_alignments"] == 1

    alignment = result["alignments"][0]

    assert alignment == {
        "aligned_sequence_1": "ATGC",
        "aligned_sequence_2": "ATGC",
        "alignment_length": 4,
        "matches": 4,
        "mismatches": 0,
        "gaps": 0,
        "gap_columns": 0,
        "identity": 1.0,
        "identity_excluding_gaps": 1.0,
        "similarity": None,
    }


def test_local_alignment_finds_best_subsequence_inside_longer_sequence():
    result = local_alignment(
        "GGATGC",
        "ATGC",
        match=2,
        mismatch=-1,
        gap_penalty=-2,
        return_all=False,
        structured=True,
    )

    assert_basic_structured_result(
        result,
        s1="GGATGC",
        s2="ATGC",
        score=8,
        best_positions=[(6, 4)],
    )

    assert result["num_alignments"] == 1

    alignment = result["alignments"][0]

    assert alignment["aligned_sequence_1"] == "ATGC"
    assert alignment["aligned_sequence_2"] == "ATGC"
    assert alignment["matches"] == 4
    assert alignment["mismatches"] == 0
    assert alignment["gaps"] == 0
    assert alignment["identity"] == 1.0


def test_local_alignment_known_biological_style_example():
    result = local_alignment(
        "GATTACA",
        "TTAC",
        match=2,
        mismatch=-1,
        gap_penalty=-2,
        return_all=False,
        structured=True,
    )

    assert_basic_structured_result(
        result,
        s1="GATTACA",
        s2="TTAC",
        score=8,
        best_positions=[(6, 4)],
    )

    assert result["num_alignments"] == 1

    alignment = result["alignments"][0]

    assert alignment["aligned_sequence_1"] == "TTAC"
    assert alignment["aligned_sequence_2"] == "TTAC"
    assert alignment["alignment_length"] == 4
    assert alignment["matches"] == 4
    assert alignment["mismatches"] == 0
    assert alignment["gaps"] == 0
    assert alignment["identity"] == 1.0


def test_local_alignment_tracks_multiple_independent_best_positions():
    result = local_alignment(
        "ATAT",
        "AT",
        match=2,
        mismatch=-1,
        gap_penalty=-2,
        return_all=True,
        structured=True,
    )

    assert_basic_structured_result(
        result,
        s1="ATAT",
        s2="AT",
        score=4,
        best_positions=[(2, 2), (4, 2)],
    )

    assert result["num_alignments"] == 2

    for alignment in result["alignments"]:
        assert alignment["aligned_sequence_1"] == "AT"
        assert alignment["aligned_sequence_2"] == "AT"
        assert alignment["matches"] == 2
        assert alignment["mismatches"] == 0
        assert alignment["gaps"] == 0
        assert alignment["identity"] == 1.0


def test_local_alignment_returns_zero_score_when_no_positive_alignment_exists():
    result = local_alignment(
        "AAAA",
        "TTTT",
        match=2,
        mismatch=-1,
        gap_penalty=-2,
        return_all=False,
        structured=True,
    )

    assert_basic_structured_result(
        result,
        s1="AAAA",
        s2="TTTT",
        score=0,
        best_positions=[],
    )

    assert result["num_alignments"] == 0
    assert result["alignments"] == []


def test_local_alignment_normalizes_positive_gap_penalty():
    result = local_alignment(
        "ATGC",
        "ATGC",
        match=2,
        mismatch=-1,
        gap_penalty=2,
        return_all=False,
        structured=True,
    )

    assert result["score"] == 8
    assert result["scoring"]["gap"] == -2
    assert result["scoring"]["gap_model"] == "linear"


def test_local_alignment_uses_custom_match_score():
    result = local_alignment(
        "AAAA",
        "AAAA",
        match=3,
        mismatch=-1,
        gap_penalty=-2,
        return_all=False,
        structured=True,
    )

    assert result["score"] == 12
    assert result["scoring"]["match"] == 3
    assert result["scoring"]["mismatch"] == -1
    assert result["scoring"]["gap"] == -2

    alignment = result["alignments"][0]

    assert alignment["aligned_sequence_1"] == "AAAA"
    assert alignment["aligned_sequence_2"] == "AAAA"
    assert alignment["matches"] == 4
    assert alignment["identity"] == 1.0


def test_local_alignment_old_return_format_still_works():
    result = local_alignment(
        "ATGC",
        "ATGC",
        match=2,
        mismatch=-1,
        gap_penalty=-2,
        return_all=False,
        structured=False,
    )

    assert result == [("ATGC", "ATGC")]

def test_local_alignment_return_all_explores_multiple_traceback_branches():
    result = local_alignment(
        "ATA",
        "ATTA",
        match=2,
        mismatch=-1,
        gap_penalty=-2,
        return_all=True,
        structured=True,
    )

    alignments = {
        (aln["aligned_sequence_1"], aln["aligned_sequence_2"])
        for aln in result["alignments"]
    }

    assert alignments == {
        ("AT", "AT"),
        ("TA", "TA"),
        ("AT-A", "ATTA"),
    }
    
def test_local_alignment_rejects_empty_sequence_input():
    with pytest.raises(ValueError):
        local_alignment("", "ATGC")

    with pytest.raises(ValueError):
        local_alignment("ATGC", "")