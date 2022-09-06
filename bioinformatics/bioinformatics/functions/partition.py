# def partition_codon_position()
# def partition_probe_flank()  # will require a function to get probe and flank (see Breinholt code)
from typing import Optional
from bioinformatics.functions.file_utils import (
    create_parent_directory,
    generate_process,
    list_files,
)


def amas_create_supermatrix(
    input_folder: str,
    partition_outfile: str,
    supermatrix_outfile: str,
    input_format: Optional[str] = "fasta",
    input_datatype: Optional[str] = "dna",
    supermatrix_format: Optional[str] = "fasta",
    partition_format: Optional[str] = "nexus",
    codons: Optional[str] = "none",
    cores: Optional[int] = 2,
) -> None:
    # https://github.com/marekborowiec/AMAS
    # supermatrix_format choices = ["fasta", "phylip", "nexus", "phylip-int", "nexus-int"]
    # partition_format choices = ["nexus", "raxml", "unspecified"]
    # codons choices = ["none", "12", "123"],
    src = "bioinformatics/src/alignment_tools/AMAS-master/amas/AMAS.py"
    # WRITE CHECK FOR input_folder CONTENT FILE TYPE
    # CALL DATATYPE CHECKER
    create_parent_directory(supermatrix_outfile)
    create_parent_directory(partition_outfile)
    cmd = []
    cmd.extend(
        [
            "python",
            src,
            "concat",
            "-d",
            input_datatype,
            "-f",
            input_format,
            "-p",
            partition_outfile,
            "-t",
            supermatrix_outfile,
            "-u",
            supermatrix_format,
            "-y",
            partition_format,
            "-n",
            codons,
            "-c",
            str(cores),
            "-i",
        ]
    )
    cmd.extend(list_files(input_folder))
    generate_process(cmd, stdout=True)
