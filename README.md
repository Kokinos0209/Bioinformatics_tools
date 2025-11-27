# Bioinformatics tools


A set of utilities for working with DNA/RNA sequences and FASTQ files.

## Functionality

**Main functions in `pybio_tools.py`:**

 `run_dna_rna_tools(*args)` function:

Performs operations on DNA/RNA sequences.

**Supported Operations:**

- `is_nucleic_acid` - validate sequence
- `transcribe` - DNA to RNA transcription  
- `reverse` - reverse sequence
- `complement` - complementary sequence
- `reverse_complement` - reverse complementary sequence

**Usage Examples:**
```python
run_dna_rna_tools('TTUU', 'is_nucleic_acid') # False 
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
```


`filter_fastq(seqs, gc_bounds, length_bounds, quality_threshold) ` function:

Filters FASTQ sequences by various criteria.

**Parameters:**

- `input_fastq - path to input FASTQ file`
- `output_fastq - path for output FASTQ file (saved in 'filtered' folder)`
- `gc_bounds - GC content bounds (default: (0, 100))`
- `length_bounds - length bounds (default: (0, 2³²))`
- `quality_threshold - average quality threshold (default: 0)`

**Usage Examples:**
```Python
filter_fastq(
    input_fastq="example_fastq.fastq",
    output_fastq="high_quality.fastq",
    quality_threshold=20,
    gc_bounds=(40, 60)
)
```

## Installation and Usage

~~~
git clone git@github.com:Kokinos0209/Bioinformatics_tools.git
cd pybio-tools
~~~

**Import functions**
~~~python
from pybio_tools import run_dna_rna_tools, filter_fastq
~~~

## Modules

**Module seq_operations**
contains functions for DNA/RNA sequence operations:
- Sequence validation
- Transcription
- Complementary transformations
- Sequence reversal

**Module fastq_operations**
contains functions for FASTQ data analysis:
- GC content calculation
- Bounds checking
- Average quality calculation
- Helper functions for filtering

## Development
This program was developed for educational purposes and can be extended with additional bioinformatics analysis functions.
Author: Elena Kokinos

