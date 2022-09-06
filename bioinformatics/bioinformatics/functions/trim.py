from typing import Optional
from clipkit import clipkit
from clipkit.modes import TrimmingMode as ClipkitMode
from bioinformatics.functions.file_utils import (
    suffix_parser,
    create_parent_directory,
    inject_grandparent_directory,
    replace_suffix,
)
from bioinformatics.models.trim import (
    TrimSoftware,
)

from bioinformatics.functions.file_utils import generate_process


def trim_alignment(
    in_file: str,
    trimmer: TrimSoftware,
    method: Optional[str] = None,
    stdout: bool = False,
) -> None:  # (construct calls to the BMGE and/or trimal programs; trimal preferred)
    cmd = []
    if method is None:
        method = "-automated1"
        out_file = inject_grandparent_directory(in_file, f"trimmed/{trimmer.name}")
    else:
        out_file = inject_grandparent_directory(
            in_file, f"trimmed/{trimmer.name}/{method}"
        )
        method = "-" + method
    create_parent_directory(out_file)
    if trimmer.name == "trimAl":
        src = "bioinformatics/src/trimming/trimal-1.4.1/source/trimal"
        in_param = "-in"
        out_param = "-out"
        # out_file = inject_prefix_suffix(out_file, "", "_trimAl")
        cmd.extend([src, in_param, in_file, out_param, out_file, method])
    elif trimmer.name == "bmge":
        src1 = "java"
        src2 = "-jar"
        src3 = "bioinformatics/src/trimming/BMGE-1.12/BMGE.jar"
        in_param = "-i"
        out_param = "-on"
        # out_file = inject_prefix_suffix(out_file, "", "_bmge")
        out_file = replace_suffix(out_file, "nex")
        type_param = "-t"
        type = "CODON"
        cmd.extend(
            [src1, src2, src3, in_param, in_file, out_param, out_file, type_param, type]
        )
    elif trimmer.name != "clipkit":
        raise Exception("Trimming software choice not understood.")
    if trimmer.name != "clipkit":
        generate_process(cmd, stdout)
    else:
        out_file = replace_suffix(out_file, "fasta")
        clipkit.execute(
            input_file=in_file,
            input_file_format=suffix_parser(in_file),
            mode=ClipkitMode("smart-gap"),
            output_file=out_file,
            output_file_format="fasta",
            gaps=0.9,
            complement=False,
            use_log=False,
        )
