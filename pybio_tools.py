from abc import ABC, abstractmethod

class BiologicalSequence(ABC):
    """Abstract base class for biological sequences"""
    def __init__(self, sequence):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)
    
    def __getitem__(self, index):
        return self.sequence[index]
    
    def __str__(self):
        return f'{self.sequence}'
 
    @abstractmethod
    def check_alphabet(self) -> bool:
        pass

class NucleicAcidSequence(BiologicalSequence):
    """Class for DNA/RNA sequences with operations like complement, reverse, and reverse complement."""

    def __init__(self, sequence):
        if type(self) == NucleicAcidSequence:
            raise NotImplementedError("Cannot instantiate NucleicAcidSequence directly")
        super().__init__(sequence)
    
    def check_alphabet(self) -> bool:
        """Check if sequence is valid DNA or RNA (contains only ATGCU and not both T and U)"""
        seq_upper = self.sequence.upper()
        right_letters = True
        for letter in seq_upper:
            if letter not in "ATGCU":
                right_letters = False
                break

        have_t = "T" in seq_upper
        have_u = "U" in seq_upper

        return right_letters and not (have_t and have_u)

    def complement(self) -> 'NucleicAcidSequence':
        """Return complementary sequence (DNA or RNA, preserving case)"""
        dna_complement_map = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
            "a": "t",
            "t": "a",
            "g": "c",
            "c": "g",
        }

        rna_complement_map = {
            "A": "U",
            "U": "A",
            "G": "C",
            "C": "G",
            "a": "u",
            "u": "a",
            "g": "c",
            "c": "g",
        }
        is_rna = "U" in self.sequence.upper()

        complement_map = rna_complement_map if is_rna else dna_complement_map
        complemented_seq = "".join([complement_map.get(nucleotide) for nucleotide in self.sequence])

        return self.__class__(complemented_seq)


    def reverse(self) -> 'NucleicAcidSequence'  :
        """Return reversed sequence"""
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self) -> 'NucleicAcidSequence':
        """Return reverse complement sequence"""
        return self.reverse().complement()


class DNASequence(NucleicAcidSequence):
    """Class for DNA sequences"""
    def __init__(self, sequence):
        super().__init__(sequence)
    
    def transcribe(self) -> 'RNASequence':
        """Transcribe DNA to RNA (replace T with U, preserving case)"""
        result = ""
        for letter in self.sequence:
            if letter == "T":
                result += "U"
            elif letter == "t":
                result += "u"
            else:
                result += letter
        return RNASequence(result)
    
class RNASequence(NucleicAcidSequence):
    """Class for RNA sequences"""
    def __init__(self, sequence):
        super().__init__(sequence)
    pass

class AminoAcidSequence(BiologicalSequence):
    """Class for amino acid sequences"""
    def __init__(self, sequence):
        super().__init__(sequence)
    
    def check_alphabet(self) -> bool:
        """Check if sequence is valid amino acid sequence (contains only 20 standard amino acids)"""
        seq_upper = self.sequence.upper()
        right_letters = True
        for letter in seq_upper:
            if letter not in "ACDEFGHIKLMNPQRSTVWY":
                right_letters = False
                break
        return right_letters
    def molecular_weight(self) -> float:
        """Calculate molecular weight (g/mol) of amino acid sequence"""
        aa_weights = {
            "A": 89.09,
            "C": 121.16,
            "D": 133.10,
            "E": 147.13,
            "F": 165.19,
            "G": 75.07,
            "H": 155.16,
            "I": 131.17,
            "K": 146.19,
            "L": 131.17,
            "M": 149.21,
            "N": 132.12,
            "P": 115.13,
            "Q": 146.15,
            "R": 174.20,
            "S": 105.09,
            "T": 119.12,
            "V": 117.15,
            "W": 204.23,
            "Y": 181.19
        }
        total_weight = sum(aa_weights.get(aa.upper(), 0) for aa in self.sequence)
        return total_weight
    

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def filter_fastq(input_fastq: str,
                 output_fastq: str,
                 gc_bounds: tuple[float, float] | float = (0, 100),
                 length_bounds: tuple[int, int] | int = (0, 2**32),
                 quality_threshold: float = 0) -> None:
    """Filters FASTQ sequences by GC content, length, and quality
    Args:
    input_fastq (str): Path to input FASTQ file
    output_fastq (str): Path to output FASTQ file 
    gc_bounds: GC content bounds (default (0, 100))
    length_bounds: Length bounds (default (0, 2**32))
    quality_threshold: Average quality threshold (default 0)

    Returns:
    dict: Filtered dictionary of sequences"""

    if isinstance(gc_bounds, (int, float)):
        gc_min, gc_max = 0.0, float(gc_bounds)
    else:
        gc_min, gc_max = gc_bounds

    if isinstance(length_bounds, int):
        len_min, len_max = 0, length_bounds
    else:
        len_min, len_max = length_bounds

    filtered_records = []
    total_records = 0

    for record in SeqIO.parse(input_fastq, "fastq"):
        total_records += 1
        gc_percent = gc_fraction(record.seq) * 100
        if not (gc_min <= gc_percent <= gc_max):
            continue

        seq_len = len(record.seq)
        if not (len_min <= seq_len <= len_max):
            continue

        qualities = record.letter_annotations["phred_quality"]
        avg_quality = sum(qualities) / len(qualities)
        if avg_quality < quality_threshold:
            continue
        filtered_records.append(record)

    print(f"Passed filter: {len(filtered_records)} reads out of {total_records}")

    SeqIO.write(filtered_records, output_fastq, "fastq")

