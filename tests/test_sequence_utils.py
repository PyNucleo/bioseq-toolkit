from bioseq.sequence_utils import (
    reverse_complement,
    gc_content,
    mrna_coding,
    mrna_template,
)


def test_reverse_complement():
    result = reverse_complement("ATGC")
    expected = "GCAT"

    assert result == expected


def test_gc_content():
    result = gc_content("GGCCAA")

    assert result == 66.66666666666666


def test_mrna_coding():
    result = mrna_coding("ATGC")

    assert result == "AUGC"


def test_mrna_template():
    result = mrna_template("ATGC")

    assert result == "UACG"