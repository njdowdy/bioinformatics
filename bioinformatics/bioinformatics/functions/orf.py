from typing import List
from Bio import SeqIO, Seq, SeqRecord
import numpy as np
from bioinformatics.functions.file_utils import (
    inject_parent_directory,
    inject_prefix_suffix,
    create_parent_directory,
)

stop_codon_code = "*"
unresolvable_codon_code = "X"

dna_code = {
    # common codes
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TGC": "C",
    "TGT": "C",
    "TGG": "W",
    # compressed notation
    "AAR": "K",
    "AAY": "N",
    "AGR": "R",
    "AGY": "S",
    "ATH": "I",
    "ATM": "I",
    "ATW": "I",
    "ATY": "I",
    "CAR": "Q",
    "CAY": "H",
    "GAR": "E",
    "GAY": "D",
    "TAY": "Y",
    "TGY": "C",
    "TTR": "L",
    "TTY": "F",
    # special compressed
    "MGR": "R",
    "RAY": "B",
    "SAR": "Z",
    "YTR": "L",
    # wobble position '--N' handled with logic
    "ACN": "T",
    "CCN": "P",
    "CGN": "R",
    "CTN": "L",
    "GCN": "A",
    "GGN": "G",
    "GTN": "V",
    "TCN": "S",
    # gap
    "---": "-",
    # stop codons
    "TRA": stop_codon_code,
    "TAR": stop_codon_code,
    "TAA": stop_codon_code,
    "TAG": stop_codon_code,
    "TGA": stop_codon_code,
}


def translate_orf(seq: Seq, orf: int, dna_code: dict[str, str]):
    if orf == 1:
        seq = "-" + seq
    elif orf == 3:
        seq = "--" + seq
    seq_triplets = [seq[i : i + 3] for i in range(0, len(seq), 3)]
    # seq_triplets_clean = [s for s in seq_triplets if len(s) == 3]
    seq_translated = [dna_code.get(u, unresolvable_codon_code) for u in seq_triplets]
    return seq_translated


def orf_score(orf_sets):
    unique_aa_per_site = [np.apply_along_axis(lambda a: len(set(a)), 0, orf_sets)][0]
    return sum(unique_aa_per_site) / len(unique_aa_per_site)


def find_best_orf(seq_list: List[str]) -> int:
    orf_scores = list(
        map(
            orf_score,
            seq_list,
        )
    )
    matches = [i for i, x in enumerate(orf_scores) if x == min(orf_scores)]
    if len(matches) == 1:
        match = matches[0] + 1
    else:
        # "Error: tied ORF choice! Selecting first occurrence (preference for forward and low valued ORFs)"
        match = matches[0] + 1
    return match


def get_best_orf(in_file: str) -> int:
    # in_file = 'bioinformatics/output/alignments/clustal_omega/LEP1/L1.fa'
    # get alignment
    aln_base = SeqIO.parse(in_file, "fasta")
    orf_list = []
    forward_orf1 = []
    reverse_orf1 = []
    forward_orf2 = []
    reverse_orf2 = []
    forward_orf3 = []
    reverse_orf3 = []
    for seq in aln_base:
        seq_data = seq.seq
        rev_seq_data = seq.seq.reverse_complement()
        forward_orf1.append(translate_orf(seq_data, 1, dna_code))
        reverse_orf1.append(translate_orf(rev_seq_data, 1, dna_code))
        forward_orf2.append(translate_orf(seq_data, 2, dna_code))
        reverse_orf2.append(translate_orf(rev_seq_data, 2, dna_code))
        forward_orf3.append(translate_orf(seq_data, 3, dna_code))
        reverse_orf3.append(translate_orf(rev_seq_data, 3, dna_code))
    orf_list.extend(
        [
            forward_orf1,
            forward_orf2,
            forward_orf3,
            reverse_orf1,
            reverse_orf2,
            reverse_orf3,
        ]
    )
    return find_best_orf(orf_list)


def fix_dna_alignment(in_file: str, best_orf: int):
    orf1_aln = []
    aln = SeqIO.parse(in_file, "fasta")
    for seq in aln:
        adj = get_adjustment_bases(seq.seq)
        orf1_aln.append(set_to_orf1(seq, best_orf, adj))
    out_file = inject_parent_directory(in_file, "ORF1")
    create_parent_directory(out_file)
    out_file = inject_prefix_suffix(out_file, "ORF1_", "")
    SeqIO.write(orf1_aln, out_file, "fasta")


def get_adjustment_bases(seq: Seq) -> tuple[str, str]:
    if seq[0] != "-":
        adj_base_start = "N"
    else:
        adj_base_start = "-"
    if seq[-1] != "-":
        adj_base_end = "N"
    else:
        adj_base_end = "-"
    return (adj_base_start, adj_base_end)


def set_to_orf1(
    seq: SeqRecord, best_orf: int, adjust_chars: tuple[str, str]
) -> SeqRecord:
    dna = seq.seq
    if best_orf in [2, 5]:  # ORF frame = 2/5, add 0 nucleotides
        dna_out = dna
    elif best_orf in [3, 6]:  # ORF frame = 3/6, add 2 undetermined nucleotides
        dna_out = adjust_chars[0] * 2 + dna
    elif best_orf in [1, 4]:  # ORF frame = 1/4, add 1 undetermined nucleotides
        dna_out = adjust_chars[0] + dna
    alnlen = len(dna_out) % 3
    if alnlen == 0:  # if length divisible by 3, do nothing
        dna_out = dna_out
    elif alnlen == 1:  # add gaps to the end such that the alignment is divisible by 3
        dna_out = dna_out + adjust_chars[1] * 2
    elif alnlen == 2:
        dna_out = dna_out + adjust_chars[1]
    seq.seq = dna_out
    return seq
