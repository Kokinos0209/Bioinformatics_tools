def is_nucleic_acid(seq: str) -> bool:
    """Check if sequence is valid DNA or RNA (contains only ATGCU and not both T and U)"""
    seq_upper = seq.upper()
    right_letters = True
    for letter in seq_upper:
        if letter not in "ATGCU":
            right_letters = False
            break

    have_t = "T" in seq_upper
    have_u = "U" in seq_upper

    return right_letters and not (have_t and have_u)
    

def transcribe(seq: str) -> str:
    """Transcribe DNA to RNA (replace T with U, preserving case)"""
    result = ""
    for letter in seq:
        if letter == "T":
            result += "U"
        elif letter == "t":
            result += "u"
        else:
            result += letter
    return result
    

def reverse(seq: str) -> str:
    """Return reversed sequence"""
    return seq[::-1]


def complement(seq: str) -> str:
    """Return complementary sequence (DNA or RNA, preserving case)"""
    dna_complement_map = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g'
    }
    
    rna_complement_map = {
        'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
        'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'
    }
    is_rna = "U" in seq.upper()

    complement_map = rna_complement_map if is_rna else dna_complement_map
    complemented_seq = ''.join([complement_map.get(nucleotide) for nucleotide in seq])
    
    return complemented_seq

def reverse_complement(seq: str) -> str:
    """Return reverse complement sequence"""
    return complement(reverse(seq))