def align_sequences(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-1):
    """Perform global alignment of two sequences using the Needleman-Wunsch algorithm."""
    n, m = len(seq1), len(seq2)
    score_matrix = [[0] * (m + 1) for _ in range(n + 1)]

    # Initialize scoring matrix
    for i in range(n + 1):
        score_matrix[i][0] = i * gap_penalty
    for j in range(m + 1):
        score_matrix[0][j] = j * gap_penalty

    # Fill the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)

    # Trace back to recover the alignment
    align1, align2 = '', ''
    i, j = n, m
    while i > 0 and j > 0:
        score_current = score_matrix[i][j]
        score_diagonal = score_matrix[i-1][j-1]
        score_left = score_matrix[i-1][j]

        if score_current == score_diagonal + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty):
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1 += seq1[i-1]
            align2 += '-'
            i -= 1
        else:
            align1 += '-'
            align2 += seq2[j-1]
            j -= 1

    # Handle trailing gaps
    while i > 0:
        align1 += seq1[i-1]
        align2 += '-'
        i -= 1
    while j > 0:
        align1 += '-'
        align2 += seq2[j-1]
        j -= 1

    return align1[::-1], align2[::-1]

def compare_sequences(ref_seq, input_seq):
    """Compare reference and input amino acid sequences and flag missense variants."""
    differences = []
    for pos, (ref_aa, input_aa) in enumerate(zip(ref_seq, input_seq), start=1):
        if ref_aa != input_aa:
            differences.append((pos, ref_aa, input_aa))
    
    return differences

def process_sequences(ref_seq, input_seq):
    """Align input sequence to reference and compare the aligned sequences."""
    aligned_ref, aligned_input = align_sequences(ref_seq, input_seq)
    return compare_sequences(aligned_ref, aligned_input)

def report_variants(differences):
    """Generate a report of variants found between aligned sequences."""
    if not differences:
        return "No variants found."
    
    report = ["Position\tReference\tInput"]
    for pos, ref, input in differences:
        report.append(f"{pos}\t{ref}\t{input}")
    
    return "\n".join(report)
