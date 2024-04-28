import pytest
from variant_tools.analysis import compare_sequences, report_missense_variants

def test_compare_sequences_identical():
    assert compare_sequences("ACDEFG", "ACDEFG") == []

def test_compare_sequences_missense():
    assert compare_sequences("ACDEFG", "ABCDEG") == [(3, 'D', 'B'), (6, 'F', 'G')]

def test_compare_sequences_error():
    with pytest.raises(ValueError):
        compare_sequences("ACDEFG", "ABCDE")

def test_report_missense_variants():
    differences = [(3, 'D', 'B'), (6, 'F', 'G')]
    expected_report = "Position\tReference\tInput\n3\tD\tB\n6\tF\tG"
    assert report_missense_variants(differences) == expected_report
