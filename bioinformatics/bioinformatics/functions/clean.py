import os
import re
from typing import Optional
from bioinformatics.functions.file_utils import inject_prefix_suffix


def standardize_data_case(
    input_file: str, case: Optional[str] = "upper", output_file: Optional[str] = None
) -> None:
    if not output_file:
        output_file = input_file + "output"  # NEEDS TO BE GENERALIZED WITH FILE UTILS
    with open(input_file) as f, open(output_file, "a") as f_out:
        for line in f:
            if line[0] != ">":
                if case == "upper":  # CHANGE TO ENUM
                    f_out.write(line.upper())
                elif case == "lower":
                    f_out.write(line.lower())
                else:
                    return print("Case option not understood!")
            else:
                f_out.write(line)


def replace_char_in_filename(
    input_file: str, search_char: str, replacement_char: Optional[str] = "_"
) -> None:
    replacement_file = input_file.replace(f"{search_char}", f"{replacement_char}")
    if replacement_file != input_file:
        os.rename(input_file, replacement_file)


def fasta_clean_newlines(input_file: str) -> None:
    output_file = input_file + ".temp"
    with open(input_file) as f, open(output_file, "w") as f_out:
        for line in f:
            if line != "\n":
                f_out.write(line)
    os.rename(output_file, input_file)


special_chars = ["[", "]", "(", ")", "{", "}", "*", "+", "?", "|", "^", "$", ".", "\\"]


def replace_ambiguous_chars(
    in_file: str, search_chars: list[str], replace_char: Optional[str] = "?"
) -> None:
    # in_file = 'bioinformatics/output/alignments/muscle3/LEP1/trimmed/bmge/ORF1_L1.nex'
    # out_file = 'bioinformatics/output/alignments/muscle3/LEP1/trimmed/bmge/ORF1_L1_fmt.nex'
    out_file = inject_prefix_suffix(in_file, "", "_fmt")
    search_chars = [f"\{x}" if x in special_chars else x for x in search_chars]
    regex = "|".join(search_chars)
    with open(in_file) as f, open(out_file, "w") as f_out:
        for line in f:
            new_line = re.sub(rf"{regex}", f"{replace_char}", line)
            f_out.write(new_line)
