import pytest
from analysis import align_sequences, compare_sequences, report_variants

def test_align_sequences():
    """Test the alignment of two sequences."""
    seq1 = "GATTACA"
    seq2 = "GCATGCU"
    aligned1, aligned2 = align_sequences(seq1, seq2)
    assert aligned1 == "G-ATTACA"
    assert aligned2 == "GCA-TGCU"

def test_compare_sequences():
    """Test comparison of aligned sequences."""
    aligned_seq1 = "G-ATTACA"
    aligned_seq2 = "GCA-TGCU"
    differences = compare_sequences(aligned_seq1, aligned_seq2)
    assert differences == [(2, '-', 'C'), (4, 'T', '-'), (6, 'A', 'G'), (7, 'C', 'C'), (8, 'A', 'U')]

def test_report_variants_no_variants():
    """Test the report output when there are no variants."""
    differences = []
    assert report_variants(differences) == "No variants found."

def test_report_variants_with_variants():
    """Test the report output when there are variants."""
    differences = [(3, 'A', 'C'), (5, 'T', 'G')]
    expected_report = "Position\tReference\tInput\n3\tA\tC\n5\tT\tG"
    assert report_variants(differences) == expected_report

