from enum import Enum
from pathlib import Path
from typing import Optional

class IndexType(Enum):
    I5 = "i5"
    I7 = "i7"

class Index:
    def __init__(self, index_type: IndexType, sequence: str):
        self.index_type = index_type
        self.sequence = sequence
        self.length = len(sequence)

class LibraryType(Enum):
    PAIRED_END = "paired-end"
    SINGLE_END = "single-end"

class SequencingLibrary:
    def __init__(
        self,
        library_type: LibraryType,
        r1_fastq: Path,
        r2_fastq: Optional[Path] = None,
        i5_index: Optional[Index] = None,
        i7_index: Optional[Index] = None
    ):
        self.library_type = library_type
        self.r1_fastq = r1_fastq
        self.r2_fastq = r2_fastq
        self.i5_index = i5_index
        self.i7_index = i7_index
        self.sequencing_stats = {}

    def to_sam(self, output_path: Path):
        pass

    def to_bam(self, output_path: Path):
        pass

    def to_cram(self, output_path: Path):
        pass

def create_seqlib(
    library_type: LibraryType,
    r1_fastq: str,
    r2_fastq: Optional[str] = None,
    i5_sequence: Optional[str] = None,
    i7_sequence: Optional[str] = None
) -> SequencingLibrary:
    r1_path = Path(r1_fastq)
    r2_path = Path(r2_fastq) if r2_fastq else None
    i5_index = Index(IndexType.I5, i5_sequence) if i5_sequence else None
    i7_index = Index(IndexType.I7, i7_sequence) if i7_sequence else None

    return SequencingLibrary(library_type, r1_path, r2_path, i5_index, i7_index)
