def gc_content(seq: str) -> float:
    """
    Calculate GC content percentage of a sequence.

    Args:
        seq: Nucleotide sequence

    Returns:
        GC content as percentage (0-100)
    """
    gc_count = seq.count("G") + seq.count("C")
    total_length = len(seq)
    return (gc_count / total_length) * 100


def check_bounds(lenght: int, bounds: int | tuple[float, float]) -> bool:
    """
    Check if value is within specified bounds.

    Args:
        value: Value to check
        bounds: Single number (upper bound) or tuple (lower, upper)

    Returns:
        True if value is within bounds, False otherwise
    """
    if type(bounds) == int:
        if lenght < bounds:
            return True
        else:
            return False
    else:
        min_lenght = bounds[0]
        max_lenght = bounds[1]
        if min_lenght <= lenght <= max_lenght:
            return True
        else:
            return False


def average_quality(q_string: str) -> float:
    """
    Calculate average quality score from Phred33 encoded string.

    Args:
        q_string: Quality string in Phred33 encoding

    Returns:
        Average quality score
    """
    numbers = []
    for symbol in q_string:
        q_score = ord(symbol) - 33
        numbers.append(q_score)

    q_total = sum(numbers)
    average = q_total / len(numbers)
    return average
