from bioseq.alignment.needleman_wunsch import global_alignment


def _reference_global_alignments(s1, s2, match, mismatch, gap_penalty):
    """Enumerate optimal alignments independently from the production traceback."""
    scores = [[0] * (len(s2) + 1) for _ in range(len(s1) + 1)]

    for row in range(1, len(s1) + 1):
        scores[row][0] = scores[row - 1][0] + gap_penalty
    for column in range(1, len(s2) + 1):
        scores[0][column] = scores[0][column - 1] + gap_penalty

    for row in range(1, len(s1) + 1):
        for column in range(1, len(s2) + 1):
            pair_score = match if s1[row - 1] == s2[column - 1] else mismatch
            scores[row][column] = max(
                scores[row - 1][column - 1] + pair_score,
                scores[row][column - 1] + gap_penalty,
                scores[row - 1][column] + gap_penalty,
            )

    def enumerate_from(row, column, aligned_1="", aligned_2=""):
        if row == 0 and column == 0:
            return [(aligned_1[::-1], aligned_2[::-1])]

        target_score = scores[row][column]
        branches = []

        if row > 0 and column > 0:
            pair_score = (
                match if s1[row - 1] == s2[column - 1] else mismatch
            )
            if scores[row - 1][column - 1] + pair_score == target_score:
                branches.extend(
                    enumerate_from(
                        row - 1,
                        column - 1,
                        aligned_1 + s1[row - 1],
                        aligned_2 + s2[column - 1],
                    )
                )

        if column > 0 and scores[row][column - 1] + gap_penalty == target_score:
            branches.extend(
                enumerate_from(
                    row,
                    column - 1,
                    aligned_1 + "-",
                    aligned_2 + s2[column - 1],
                )
            )

        if row > 0 and scores[row - 1][column] + gap_penalty == target_score:
            branches.extend(
                enumerate_from(
                    row - 1,
                    column,
                    aligned_1 + s1[row - 1],
                    aligned_2 + "-",
                )
            )

        return branches

    return scores[-1][-1], enumerate_from(len(s1), len(s2))


def test_global_alignment_return_all_retains_each_tied_branch_in_order():
    expected = [
        ("-A", "AA"),
        ("A-", "AA"),
    ]

    first_result = global_alignment(
        "A",
        "AA",
        match=1,
        mismatch=-1,
        gap_penalty=-1,
        return_all=True,
        structured=False,
    )
    second_result = global_alignment(
        "A",
        "AA",
        match=1,
        mismatch=-1,
        gap_penalty=-1,
        return_all=True,
        structured=False,
    )

    assert first_result == expected
    assert second_result == expected


def test_global_alignment_return_all_retains_reverse_orientation_ties_in_order():
    result = global_alignment(
        "AA",
        "A",
        match=1,
        mismatch=-1,
        gap_penalty=-1,
        return_all=True,
        structured=False,
    )

    assert result == [
        ("AA", "-A"),
        ("AA", "A-"),
    ]


def test_global_alignment_single_traceback_keeps_first_tied_movement():
    result = global_alignment(
        "A",
        "AA",
        match=1,
        mismatch=-1,
        gap_penalty=-1,
        return_all=False,
        structured=False,
    )

    assert result == [("-A", "AA")]


def test_global_alignment_return_all_structured_tied_results_are_consistent():
    result = global_alignment(
        "A",
        "AA",
        match=1,
        mismatch=-1,
        gap_penalty=-1,
        return_all=True,
        structured=True,
    )

    assert result["algorithm"] == "Needleman-Wunsch"
    assert result["mode"] == "global"
    assert result["sequence_1"] == "A"
    assert result["sequence_2"] == "AA"
    assert result["score"] == 0
    assert result["scoring"] == {
        "match": 1,
        "mismatch": -1,
        "gap_penalty": -1,
        "matrix": None,
        "gap_model": "linear",
    }
    assert result["num_alignments"] == 2
    assert result["num_alignments"] == len(result["alignments"])
    assert result["alignments"] == [
        {
            "aligned_sequence_1": "-A",
            "aligned_sequence_2": "AA",
            "alignment_length": 2,
            "matches": 1,
            "mismatches": 0,
            "gaps": 1,
            "gap_columns": 1,
            "identity": 0.5,
            "identity_excluding_gaps": 1.0,
            "similarity": None,
        },
        {
            "aligned_sequence_1": "A-",
            "aligned_sequence_2": "AA",
            "alignment_length": 2,
            "matches": 1,
            "mismatches": 0,
            "gaps": 1,
            "gap_columns": 1,
            "identity": 0.5,
            "identity_excluding_gaps": 1.0,
            "similarity": None,
        },
    ]


def test_global_alignment_return_all_accumulates_nested_tie_leaves():
    expected = [
        ("--A", "AAA"),
        ("-A-", "AAA"),
        ("A--", "AAA"),
    ]

    # The optimal path branches at (1, 3); its left branch reaches another tie
    # at (1, 2), so retaining every complete leaf requires nested accumulation.
    reference_score, reference_alignments = _reference_global_alignments(
        "A", "AAA", match=1, mismatch=-1, gap_penalty=-1
    )
    result = global_alignment(
        "A",
        "AAA",
        match=1,
        mismatch=-1,
        gap_penalty=-1,
        return_all=True,
        structured=False,
    )

    assert reference_score == -1
    assert reference_alignments == expected
    assert result == expected


def test_global_alignment_no_tie_control_has_one_alignment_in_each_format():
    all_alignments = global_alignment(
        "AC",
        "AC",
        match=1,
        mismatch=-1,
        gap_penalty=-1,
        return_all=True,
        structured=False,
    )
    single_alignment = global_alignment(
        "AC",
        "AC",
        match=1,
        mismatch=-1,
        gap_penalty=-1,
        return_all=False,
        structured=False,
    )
    structured_result = global_alignment(
        "AC",
        "AC",
        match=1,
        mismatch=-1,
        gap_penalty=-1,
        return_all=True,
        structured=True,
    )

    assert all_alignments == [("AC", "AC")]
    assert single_alignment == all_alignments
    assert structured_result["num_alignments"] == 1
    assert len(structured_result["alignments"]) == 1
