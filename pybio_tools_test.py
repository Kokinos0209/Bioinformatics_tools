import pytest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pybio_tools import filter_fastq


def test_gc_filter_keeps_right_gc_read(tmp_path):
    input_file = tmp_path / "input.fastq"
    output_file = tmp_path / "output.fastq"

    record1 = SeqRecord(Seq("GCGCGCA"), id="good_gc", description="tessting read1")
    record1.letter_annotations["phred_quality"] = [30, 30, 30, 30, 30, 30, 30]

    record2 = SeqRecord(Seq("ATATG"), id="bad_gc", description="tessting read2")
    record2.letter_annotations["phred_quality"] = [30, 30, 30, 30, 30]
    SeqIO.write([record1, record2], input_file, "fastq")

    filter_fastq(str(input_file), str(output_file), gc_bounds=(80, 100))

    result = list(SeqIO.parse(output_file, "fastq"))

    assert result[0].id == "good_gc"


def test_length_filter_keeps_right_length(tmp_path):
    input_file = tmp_path / "input.fastq"
    output_file = tmp_path / "output.fastq"

    record1 = SeqRecord(Seq("GCGCGCA"), id="len7", description="tessting read1")
    record1.letter_annotations["phred_quality"] = [30, 30, 30, 30, 30, 30, 30]

    record2 = SeqRecord(Seq("ATATG"), id="len5", description="tessting read2")
    record2.letter_annotations["phred_quality"] = [30, 30, 30, 30, 30]

    SeqIO.write([record1, record2], input_file, "fastq")

    filter_fastq(str(input_file), str(output_file), length_bounds=(5, 6))

    result = list(SeqIO.parse(output_file, "fastq"))

    assert result[0].id == "len5"



def test_quality_filter_keeps_good_quality(tmp_path):
    input_file = tmp_path / "input.fastq"
    output_file = tmp_path / "output.fastq"

    record1 = SeqRecord(Seq("ATGC"), id="good_q", description="")
    record1.letter_annotations["phred_quality"] = [40, 40, 40, 40]

    record2 = SeqRecord(Seq("ATGC"), id="bad_q", description="")
    record2.letter_annotations["phred_quality"] = [10, 10, 10, 10]

    SeqIO.write([record1, record2], input_file, "fastq")

    filter_fastq(str(input_file), str(output_file), quality_threshold=30)

    result = list(SeqIO.parse(output_file, "fastq"))

    assert result[0].id == "good_q"


def test_quality_filter_removes_bad_quality(tmp_path):
    input_file = tmp_path / "input.fastq"
    output_file = tmp_path / "output.fastq"

    record = SeqRecord(Seq("ATGC"), id="bad_q", description="")
    record.letter_annotations["phred_quality"] = [10, 10, 10, 10]

    SeqIO.write([record], input_file, "fastq")

    filter_fastq(str(input_file), str(output_file), quality_threshold=20)

    result = list(SeqIO.parse(output_file, "fastq"))

    assert len(result) == 0


def test_output_file_is_created(tmp_path):
    input_file = tmp_path / "input.fastq"
    output_file = tmp_path / "output.fastq"

    record = SeqRecord(Seq("ATGC"), id="read1", description="")
    record.letter_annotations["phred_quality"] = [30, 30, 30, 30]

    SeqIO.write([record], input_file, "fastq")

    filter_fastq(str(input_file), str(output_file))

    assert output_file.exists()


def test_error_if_input_file_does_not_exist(tmp_path):
    input_file = tmp_path / "missing.fastq"
    output_file = tmp_path / "output.fastq"

    with pytest.raises(FileNotFoundError):
        filter_fastq(str(input_file), str(output_file))


