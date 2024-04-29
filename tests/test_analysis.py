import pytest
from variant_tools.analysis import compare_sequences, report_missense_variants

def test_compare_sequences_identical():
    """Test that no differences are found when sequences are identical."""
    assert compare_sequences("ACDEFG", "ACDEFG") == []

def test_compare_sequences_missense():
    """Test that correct differences are identified."""
    assert compare_sequences("ACDEFG", "ABCDEG") == [(3, 'D', 'B'), (6, 'F', 'G')]

def test_compare_sequences_error():
    """Test that an error is raised when sequences are of different lengths."""
    with pytest.raises(ValueError):
        compare_sequences("ACDEFG", "ABCDE")

def test_report_missense_variants_no_variants():
    """Test the report output when there are no variants."""
    differences = []
    assert report_missense_variants(differences) == "No missense variants found."

def test_report_missense_variants_with_variants():
    """Test the report output when there are variants."""
    differences = [(3, 'D', 'B'), (6, 'F', 'G')]
    expected_report = "Position\tReference\tInput\n3\tD\tB\n6\tF\tG"
    assert report_missense_variants(differences) == expected_report

