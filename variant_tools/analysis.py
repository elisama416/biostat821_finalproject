def compare_sequences(ref_seq, input_seq):
    """Compare reference and input amino acid sequences and flag missense variants."""
    if len(ref_seq) != len(input_seq):
        raise ValueError("Sequences must be of the same length.")

    differences = []
    for pos, (ref_aa, input_aa) in enumerate(zip(ref_seq, input_seq), start=1):
        if ref_aa != input_aa:
            differences.append((pos, ref_aa, input_aa))
    
    return differences

def report_missense_variants(differences):
    """Generate a report of missense variants."""
    if not differences:
        return "No missense variants found."
    
    report = ["Position\tReference\tInput"]
    for pos, ref, input in differences:
        report.append(f"{pos}\t{ref}\t{input}")
    
    return "\n".join(report)
