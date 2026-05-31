from bioseq.alignment.needleman_wunsch import global_alignment


def test_global_alignment_structured_result_with_default_gap_penalty():
    result = global_alignment(
        "ATGCG",
        "ATCGA",
        match=1,
        mismatch=-1,
        gap_penalty=-2,
        return_all=False,
        structured=True,
    )

    assert result["algorithm"] == "Needleman-Wunsch"
    assert result["mode"] == "global"
    assert result["sequence_1"] == "ATGCG"
    assert result["sequence_2"] == "ATCGA"
    assert result["score"] == 0

    assert result["scoring"] == {
        "match": 1,
        "mismatch": -1,
        "gap": -2,
        "matrix": None,
        "gap_model": "linear",
    }

    assert result["num_alignments"] == 1
    assert isinstance(result["alignments"], list)

    alignment = result["alignments"][0]

    assert alignment == {
        "aligned_sequence_1": "ATGCG-",
        "aligned_sequence_2": "AT-CGA",
        "alignment_length": 6,
        "matches": 4,
        "mismatches": 0,
        "gaps": 2,
        "gap_columns": 2,
        "identity": 0.667,
        "identity_excluding_gaps": 1.0,
        "similarity": None,
    }


def test_global_alignment_structured_result_with_lower_gap_penalty():
    result = global_alignment(
        "ATGCG",
        "ATCGA",
        match=1,
        mismatch=-1,
        gap_penalty=-1,
        return_all=False,
        structured=True,
    )

    assert result["score"] == 2
    assert result["scoring"]["gap"] == -1

    alignment = result["alignments"][0]

    assert alignment["aligned_sequence_1"] == "ATGCG-"
    assert alignment["aligned_sequence_2"] == "AT-CGA"
    assert alignment["matches"] == 4
    assert alignment["mismatches"] == 0
    assert alignment["gaps"] == 2
    assert alignment["identity"] == 0.667
    assert alignment["identity_excluding_gaps"] == 1.0


def test_global_alignment_structured_result_with_higher_match_score():
    result = global_alignment(
        "ATGCG",
        "ATCGA",
        match=2,
        mismatch=-1,
        gap_penalty=-1,
        return_all=False,
        structured=True,
    )

    assert result["score"] == 6
    assert result["scoring"]["match"] == 2
    assert result["scoring"]["mismatch"] == -1
    assert result["scoring"]["gap"] == -1

    alignment = result["alignments"][0]

    assert alignment["aligned_sequence_1"] == "ATGCG-"
    assert alignment["aligned_sequence_2"] == "AT-CGA"
    assert alignment["alignment_length"] == 6
    assert alignment["matches"] == 4
    assert alignment["mismatches"] == 0
    assert alignment["gaps"] == 2
    assert alignment["gap_columns"] == 2


def test_global_alignment_structured_result_for_identical_sequences():
    result = global_alignment(
        "GATTACA",
        "GATTACA",
        match=1,
        mismatch=-1,
        gap_penalty=-2,
        return_all=False,
        structured=True,
    )

    assert result["score"] == 7
    assert result["num_alignments"] == 1

    alignment = result["alignments"][0]

    assert alignment["aligned_sequence_1"] == "GATTACA"
    assert alignment["aligned_sequence_2"] == "GATTACA"
    assert alignment["alignment_length"] == 7
    assert alignment["matches"] == 7
    assert alignment["mismatches"] == 0
    assert alignment["gaps"] == 0
    assert alignment["gap_columns"] == 0
    assert alignment["identity"] == 1.0
    assert alignment["identity_excluding_gaps"] == 1.0


def test_global_alignment_old_return_format_still_works():
    result = global_alignment(
        "GATTACA",
        "GATTACA",
        return_all=False,
        structured=False
    )

    assert result == [("GATTACA", "GATTACA")]