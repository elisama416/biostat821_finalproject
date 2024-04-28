import argparse
from variant_tools.analysis import compare_sequences, report_missense_variants

def main():
    parser = argparse.ArgumentParser(description="Compare amino acid sequences for missense mutations.")
    parser.add_argument('ref_seq', type=str, help='Reference amino acid sequence.')
    parser.add_argument('input_seq', type=str, help='Input amino acid sequence to compare.')
    args = parser.parse_args()

    try:
        differences = compare_sequences(args.ref_seq, args.input_seq)
        report = report_missense_variants(differences)
        print(report)
    except ValueError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
