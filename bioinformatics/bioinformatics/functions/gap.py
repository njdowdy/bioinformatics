import re
from Bio import SeqIO
from bioinformatics.functions.file_utils import (
    suffix_parser,
)


def compile_gap_pattern(gap_chars: list[str]):
    pattern = ""
    for char in gap_chars:
        pattern = pattern + f"{char}+|"
    pattern = re.sub(r"\|$", "", pattern)
    pattern = pattern.replace("?", "\?")
    return pattern


def find_unique_gaps(
    in_file: str,
    gap_chars: list[str] = ["-", "X", "?"],
    score_together: bool = True,
    filter: list[str] = ["start", "end", "center"],
    ignore_trailing_ambigs: bool = True,
) -> set[tuple[int, int]]:
    # in_file = 'bioinformatics/output/alignments/clustal_omega/LEP1/L1.fa'
    # in_file = 'bioinformatics/output/alignments/muscle3/LEP1/trimmed/augment/ORF1_L1_fmt.nex'
    if suffix_parser(in_file) == "phylip":
        file_format = "phylip-relaxed"
    else:
        file_format = suffix_parser(in_file)
    aln_base = SeqIO.parse(in_file, file_format)
    gaps = set()
    if score_together:
        pattern = "|".join(gap_chars).replace("?", "\?")
    else:
        pattern = compile_gap_pattern(gap_chars)
    for seq in aln_base:
        sequence = str(seq.seq)
        if ignore_trailing_ambigs:
            subpattern = ("|".join(gap_chars) + "{1,3}").replace("?", "\\?")
            sequence = re.sub(rf"{subpattern}$", "", sequence)
        gap_list = [x.span() for x in re.finditer(rf"[{pattern}]+", sequence)]
        if gap_list:
            gaps = gaps | set(gap_list)
    gaps_out = set()
    if "start" in filter:
        gaps_out = gaps_out | {s for s in list(gaps) if s[0] == 0}
    if "center" in filter:
        gaps_out = gaps_out | {
            s for s in list(gaps) if s[0] != 0 and s[1] != len(seq.seq)
        }
    if "end" in filter:
        gaps_out = gaps_out | {s for s in list(gaps) if s[1] == len(seq.seq)}
    gaps_out = sorted(gaps_out)
    return gaps_out


def score_gaps(in_file: str, gaps: set[tuple[int, int]]):
    score_matrix = {}
    if suffix_parser(in_file) == "phylip":
        file_format = "phylip-relaxed"
    else:
        file_format = suffix_parser(in_file)
    aln_base = SeqIO.parse(in_file, file_format)
    # for gap in gaps:
    #     for seq in aln_base:
    # check if gap overlaps the next gap
    return score_matrix
