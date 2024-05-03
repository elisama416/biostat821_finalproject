# Genetic Sequence Analysis Library

The Genetic Sequence Analysis Library is designed to assist with the alignment, comparison, and analysis of genetic sequences, particularly focused on identifying missense variants from amino acid sequences.

## For End Users

### Setup/Installation Instructions

To get started with the Genetic Sequence Analysis Library, follow these steps:

1. Clone the repository:
   ```python
   git clone https://github.com/elisama416/biostat821_finalproject.git
   ```
2. Install the required dependencies
   ```python
   pip install -r requirements-test.txt
   ```

The library expects sequences to be provided as strings in the amino acid single-letter format.

### Examples

Here are some examples of how to use the library:
```python
from variant_tools.analysis import align_sequences, compare_sequences, process_sequences, report_variants

# Align two sequences
aligned_seq1, aligned_seq2 = align_sequences("GATTACA", "GCATGCU")

# Compare aligned sequences
differences = compare_sequences(aligned_seq1, aligned_seq2)

# Generate a report of variants found
variant_report = report_variants(differences)
print(variant_report)
```

## For Contributors

To contribute to the EHR Library, ensure that you have 'pytest' installed:
```python
pip install pytest
````

Run the tests locally from the project root directory:
```python
pytest
```

Ensure that all tests pass before submitting a pull request. Here's an example of how a test might look for the 'parse_date' function:
```python
import pytest
from variant_tools.analysis import align_sequences

def test_align_sequences_basic():
    seq1 = "GATTACA"
    seq2 = "GCATGCU"
    aligned1, aligned2 = align_sequences(seq1, seq2)
    assert aligned1 == "G-ATTACA" and aligned2 == "GCA-TGCU"
```

Refer to the 'tests/' directory for more examples.

