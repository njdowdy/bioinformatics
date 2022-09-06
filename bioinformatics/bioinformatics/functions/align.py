from typing import Optional
from Bio import AlignIO, Seq
from bioinformatics.functions.file_utils import (
    replace_parent_directory,
    suffix_parser,
    create_parent_directory,
    inject_parent_directory,
    replace_suffix,
)
from bioinformatics.models.alignment import (
    AlignmentOutputFormat,
    AlignmentSoftware,
)
from bioinformatics.functions.file_utils import generate_process
from bioinformatics.functions.ingest import read_seq_file, write_seq_file


def pad_alignment(in_file: str, ambiguious_char: Optional[str] = "-") -> None:
    suffix = suffix_parser(in_file)
    records = list(read_seq_file(in_file))
    maxlen = max(len(record.seq) for record in records)
    for record in records:
        if len(record.seq) != maxlen:
            sequence = str(record.seq).ljust(maxlen, ambiguious_char)
            record.seq = Seq.Seq(sequence)
    output_file = inject_parent_directory(in_file, "padded")
    output_file = replace_suffix(output_file, suffix)
    create_parent_directory(output_file)
    write_seq_file(records, output_file, suffix)


def align_pairwise(
    fasta: str,
):
    return None


def format_parser(in_file_suffix: str) -> str:
    if in_file_suffix in ["fa", "fas", "fasta"]:
        in_file_format = "fasta"
    elif in_file_suffix in ["phy", "phylip"]:
        in_file_format = "phylip"
    else:
        raise Exception("File format could not be parsed. Please provide.")
    return in_file_format


def convert_format(
    in_file: str,
    out_file: str,
    in_file_format: Optional[str] = None,
    out_file_format: Optional[str] = None,
) -> None:
    in_file_suffix = suffix_parser(in_file)
    if not in_file_format:
        in_file_format = format_parser(in_file_suffix)
    out_file_suffix = suffix_parser(out_file)
    if not out_file_format:
        out_file_format = format_parser(out_file_suffix)
    create_parent_directory(out_file)
    alignments = AlignIO.parse(in_file, in_file_format)
    AlignIO.write(alignments, open(out_file, "w"), out_file_format)


def perform_alignment(
    in_file: str,
    aligner: AlignmentSoftware,
    output_format: AlignmentOutputFormat,
    iterations: Optional[int] = None,
    stdout: bool = False,
) -> None:
    # in_file = 'bioinformatics/input/fasta/LEP1/L1.fa'
    # aligner = AlignmentSoftware("clustal_omega")
    # output_format = AlignmentOutputFormat("fasta")
    # iterations = 5

    output_format = output_format.name
    out_file = replace_parent_directory(
        in_file, "/".join(in_file.split("/")[1:3]), f"output/alignments/{aligner.name}"
    )
    out_file = out_file.replace(suffix_parser(out_file), output_format)
    create_parent_directory(out_file)
    if not iterations:
        iterations = 1
    cmd = []
    if aligner.name == "clustal_omega":
        src = "bioinformatics/src/align/clustalo-1.2.4-Ubuntu-x86_64"
        in_param = "-i"
        out_param = "-o"
        iters_param = "--iter"
        cmd.extend(
            [
                src,
                in_param,
                in_file,
                out_param,
                out_file,
                "--outfmt",
                output_format,
                iters_param,
                str(iterations),
            ]
        )
    elif aligner.name == "muscle5":
        src = "bioinformatics/src/align/muscle5.1.linux_intel64"
        in_param = "-align"
        out_param = "-output"
        iters_param = "-replicates"
        cmd.extend(
            [src, in_param, in_file, out_param, out_file, iters_param, str(iterations)]
        )
    elif aligner.name == "muscle3":
        src = "bioinformatics/src/align/muscle3.8.31_i86linux64"
        in_param = "-in"
        out_param = "-out"
        iters_param = "-maxiters"
        cmd.extend(
            [src, in_param, in_file, out_param, out_file, iters_param, str(iterations)]
        )
    else:
        raise Exception("Alignment software choice not understood.")
    generate_process(cmd, stdout)
