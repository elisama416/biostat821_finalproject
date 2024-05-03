"""This tests the functions from the variant analysis toolset."""

from typing import List, Tuple

from variant_tools.analysis import (
    align_sequences,
    compare_sequences,
    process_sequences,
    report_variants,
)


def test_align_sequences() -> None:
    """Test the alignment of two sequences."""
    seq1 = "GATTACA"
    seq2 = "GCATGCU"
    aligned1, aligned2 = align_sequences(seq1, seq2)
    assert aligned1 == "G-ATTACA", f"Expected 'G-ATTACA', got '{aligned1}'"
    assert aligned2 == "GCA-TGCU", f"Expected 'GCA-TGCU', got '{aligned2}'"

def test_compare_sequences() -> None:
    """Test comparison of aligned sequences."""
    aligned_seq1 = "G-ATTACA"
    aligned_seq2 = "GCA-TGCU"
    differences = compare_sequences(aligned_seq1, aligned_seq2)
    expected_differences = [(2, '-', 'C'), (4, 'T', '-'), (6, 'A', 'G'), \
                            (8, 'A', 'U')]
    assert differences == expected_differences, f"Expected \
        {expected_differences}, got {differences}"

def test_process_sequences() -> None:
    """Test the full process from alignment to comparison."""
    ref_seq = "GATTACA"
    input_seq = "GCATGCU"
    differences = process_sequences(ref_seq, input_seq)
    expected_differences = [(2, '-', 'C'), (4, 'T', '-'), (6, 'A', 'G'), \
                            (8, 'A', 'U')]
    assert differences == expected_differences, f"Expected \
        {expected_differences}, got {differences}"

def test_report_variants_no_variants() -> None:
    """Test the report output when there are no variants."""
    differences: List[Tuple[int, str, str]] = []
    assert report_variants(differences) == "No variants found.", "Report \
        should indicate no variants found"

def test_report_variants_with_variants() -> None:
    """Test the report output when there are variants."""
    differences = [(3, 'A', 'C'), (5, 'T', 'G')]
    expected_report = "Position\tReference\tInput\n3\tA\tC\n5\tT\tG"
    assert report_variants(differences) == expected_report, f"Expected \
        '{expected_report}', got '{report_variants(differences)}'"

def test_align_sequences_different_lengths() -> None:
    """Test the alignment of sequences of different lengths."""
    seq1 = "GATTACA"
    seq2 = "GCAT"
    aligned1, aligned2 = align_sequences(seq1, seq2)
    expected_aligned1 = "GATTACA-"
    expected_aligned2 = "----GCAT"
    assert aligned1 == expected_aligned1, f"Expected '{expected_aligned1}', \
        got '{aligned1}'"
    assert aligned2 == expected_aligned2, f"Expected '{expected_aligned2}', \
        got '{aligned2}'"


