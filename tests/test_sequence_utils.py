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

    assert result == 66.67


def test_mrna_coding():
    result = mrna_coding("ATGC")

    assert result == "AUGC"


def test_mrna_template():
    result = mrna_template("ATGC")

    assert result == "UACG"

def test_sequence_utils_accept_lowercase_dna():
    assert reverse_complement("atgc") == "GCAT"
    assert gc_content("ggccaa") == 66.67
    assert mrna_coding("atgc") == "AUGC"
    assert mrna_template("atgc") == "UACG"