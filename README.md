# Bioinformatics tools

![alt text](image.png)

## Features

- Object-oriented API for biological sequences:
	- DNASequence/RNASequence with complement, reverse and reverse-complement operations
	- Transcription of DNA to RNA
	- AminoAcidSequence with molecular weight calculation
- FASTQ filtering by GC content, read length and average quality (using Biopython).

## Installation

Make sure you have Python 3 and Biopython installed:

```bash
pip install biopython
```

Clone the repository and move into the project directory:

```bash
git clone git@github.com:Kokinos0209/Bioinformatics_tools.git
cd pybio-tools
```

## Usage examples

### Working with DNA and RNA

```python
from pybio_tools import DNASequence, RNASequence

seq_dna = DNASequence("ATGCGTAG")
print(seq_dna)                 # ATGCGTAG
seq_dna.check_alphabet()       # True

seq_rna = seq_dna.transcribe() # DNA -> RNA
print(seq_rna)                 # AUGCGUAG

seq_dna_complement = seq_dna.complement() #TACGCATC
seq_dna_reverse = seq_dna.reverse() #GATGCGTA
seq_dna_rev_comp = seq_dna.reverse_complement() #CTACGCAT
```
P.S. always make check_alphabet for your seq!

### Amino acid sequences

```python
from pybio_tools import AminoAcidSequence

seq_AA = AminoAcidSequence("ACDEFGHIKLMNPQRSTVWY")
seq_AA.check_alphabet()            # True
seq_AA_weight = seq_AA.molecular_weight()
print(f"{seq_AA_weight:.2f} g/mol") # 2738.02 g/mol   
```
P.S. always make check_alphabet for your seq!

### FASTQ filtering

```python
from pybio_tools import filter_fastq

filter_fastq(
		input_fastq="input.fastq",
		output_fastq="filtered.fastq",
		gc_bounds=(40, 60),        # GC% range
		length_bounds=(50, 300),   # read length range
		quality_threshold=20       # minimum average quality
)
```

The function prints how many reads passed the filter and writes them to the output FASTQ file.

## Development
This program was developed for educational purposes and can be extended with additional bioinformatics analysis functions.
Author: Elena Kokinos

