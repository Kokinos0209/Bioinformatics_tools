from modules import seq_operations, fastq_operations


def run_dna_rna_tools(*args):
    """
    Perform DNA/RNA sequence operations

    Args:
        *args: Sequences (str) and operation (last arg).
               Operations: 'is_nucleic_acid', 'transcribe', 'reverse',
                           'complement', 'reverse_complement'

    Returns:
        Single result if one sequence, list if multiple
        None for invalid sequences
        Preserves letter case
    """
    *seqs, operation = args
    results = []
    for seq in seqs:
        if operation == "is_nucleic_acid":
            results.append(seq_operations.is_nucleic_acid(seq))
        else:
            if not seq_operations.is_nucleic_acid(seq):
                results.append(None)
            else:
                if operation == "transcribe":
                    results.append(seq_operations.transcribe(seq))
                elif operation == "reverse":
                    results.append(seq_operations.reverse(seq))
                elif operation == "complement":
                    results.append(seq_operations.complement(seq))
                elif operation == "reverse_complement":
                    results.append(seq_operations.reverse_complement(seq))

    return results


def filter_fastq(
    seqs: dict[str, tuple[str, str]],
    gc_bounds: tuple[float, float] | float = (0, 100),
    length_bounds: tuple[int, int] | int = (0, 2**32),
    quality_threshold: float = 0,
) -> dict[str, tuple[str, str]]:
    """
    Filters FASTQ sequences by GC content, length, and quality

    Args:
    seqs (dict): A dictionary {sequence_name: (sequence, quality)}
    gc_bounds: GC content bounds (default (0, 100))
    length_bounds: Length bounds (default (0, 2**32))
    quality_threshold: Average quality threshold (default 0)

    Returns:
    dict: Filtered dictionary of sequences
    """
    result = {}
    for name in seqs:
        sequence_data = seqs[name]
        sequence = sequence_data[0]
        quality = sequence_data[1]

        gc_percent = fastq_operations.gc_content(sequence)
        gc_opt = fastq_operations.check_bounds(gc_percent, gc_bounds)
        if not gc_opt:
            continue

        length_seq = len(sequence)
        length_opt = fastq_operations.check_bounds(length_seq, length_bounds)
        if not length_opt:
            continue

        avg_quality = fastq_operations.average_quality(quality)
        if avg_quality < quality_threshold:
            continue

        result[name] = (sequence, quality)

    return result
