from io import TextIOWrapper
from typing import Optional, List
from Bio import SeqIO, SeqRecord
from Bio.SeqIO.FastaIO import FastaIterator
from bioinformatics.functions.file_utils import (
    parse_output_file_name,
    suffix_parser,
)


def read_seq_file(in_file: str) -> SeqIO:
    suffix = suffix_parser(in_file)
    records = SeqIO.parse(in_file, suffix)
    return records


def write_seq_file(
    records: list[SeqRecord.SeqRecord], out_file: str, suffix: str
) -> None:
    SeqIO.write(records, out_file, suffix)


def substring_label_search(
    record: FastaIterator, filter: List[str], output_file: TextIOWrapper
) -> bool:
    discovered = False
    for substring in filter:
        if substring in record.description:
            output_file.write(f">{record.id}\n{record.seq}\n")
            discovered = True
    return discovered


def filter_fasta_by_label(
    fasta: str, primary_filter: List[str], secondary_filter: Optional[List[str]]
):
    output_file = parse_output_file_name(fasta, "output/taxon_filtered_alignments")
    lab_found = False
    with open(output_file, "w") as f:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            discovered = substring_label_search(seq_record, primary_filter, f)
            if discovered:
                lab_found = True
        if not lab_found:
            for seq_record in SeqIO.parse(fasta, "fasta"):
                discovered = substring_label_search(seq_record, secondary_filter, f)
                if discovered:
                    lab_found = True
        #     if lab_found:
        #         print(f"Substitute Was Required & Found For: {output_file}")
        #     else:
        #         print(f"No Substitute Available For: {output_file}")
        # else:
        #     print(f"Taxa Grabbed for: {output_file}")


def get_all_tips(fasta: str, tip_list: str = []):
    for seq_record in SeqIO.parse(fasta, "fasta"):
        locus_name = ".".join(fasta.split("/")[-1].split(".")[0:-1])
        name = seq_record.description.replace(f"{locus_name}", "")
        if not name in tip_list:
            tip_list.append(name)
    return tip_list


def get_molecular_data_type(input_file: str) -> str:
    nucleotides = ["A", "C", "T", "G"]
    amino_acids = [
        "R",
        "D",
        "C",
        "E",
        "Q",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
    ]
    with open(input_file) as f:
        for line in f:
            try:
                if ">" not in line:
                    if any(nucleotide in line.upper() for nucleotide in nucleotides):
                        return "nucl"
                    elif any(amino_acid in line.upper() for amino_acid in amino_acids):
                        return "Protein"
                    else:
                        raise ValueError("Could not parse data type")
            except ValueError:
                raise
