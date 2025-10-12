import os

def read_fastq(input_file):
    sequences = {}
    with open(input_file, 'r') as file:
        lines = file.readlines()
    
    for i in range(0, len(lines), 4):
        header = lines[i].strip()
        sequence = lines[i+1].strip()
        plus = lines[i+2].strip()
        quality = lines[i+3].strip()
        
        name = header[1:]
        sequences[name] = (sequence, quality)
    
    return sequences

def write_fastq(sequences_dict, output_file):
    folder = os.path.dirname(output_file)
    if folder and not os.path.exists(folder):
        os.makedirs(folder)
        
    if os.path.exists(output_file):
        raise FileExistsError(f"File {output_file} already exists. Will not overwrite)))")
    
    with open(output_file, 'w') as file:
        for name, (sequence, quality) in sequences_dict.items():
            file.write(f"@{name}\n")
            file.write(f"{sequence}\n")
            file.write("+\n")
            file.write(f"{quality}\n")

def check_gc_content(seq: str) -> float:
    """
    Calculate GC content percentage of a sequence.

    Args:
        seq: Nucleotide sequence

    Returns:
        GC content as percentage (0-100)
    """
    seq_upper = seq.upper()
    gc_count = seq_upper.count("G") + seq_upper.count("C")
    total_length = len(seq_upper)
    return (gc_count / total_length) * 100


def check_bounds(value: int, bounds: int | tuple[float, float]) -> bool:
    """
    Check if value is within specified bounds.

    Args:
        value: Value to check
        bounds: Single number (upper bound) or tuple (lower, upper)

    Returns:
        True if value is within bounds, False otherwise
    """
    if type(bounds) == int:
         return value < bounds
    else:
        min_value = bounds[0]
        max_value = bounds[1]
        return min_value <= value <= max_value


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
