from io import TextIOWrapper
import re
from typing import Optional, List, TextIO
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.SeqIO import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo
from bioinformatics.functions.file_utils import create_parent_directory, suffix_parser


def parse_output_file_name(in_file: str, output_folder: str = "output") -> str:
    suffix = suffix_parser(in_file)
    output_file = re.sub(rf"{suffix}$", f".out{suffix}", in_file).replace("..", ".")
    output_file = re.sub(r"/input/fasta/", f"/{output_folder}/", output_file)
    create_parent_directory(output_file)
    return output_file


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


def generate_consensus(fasta: str, threshold: Optional[float] = 0.5):
    # fasta = 'bioinformatics/output/taxon_filtered_alignments/ZFMK1/EOG090R0A51_3.final.out.fas'
    output_file = parse_output_file_name(fasta).replace(".out.out", ".out")
    subfolders_to_replace = "/".join(
        output_file.split("/")[output_file.split("/").index("output") + 1 : -2]
    )
    output_file = output_file.replace(f"{subfolders_to_replace}", "consensus_sequences")
    locus_name = (
        re.search(r"consensus_sequences/(.*?)\.out", output_file)
        .group(1)
        .replace("/", ".")
    )
    try:
        alignment = AlignIO.read(fasta, "fasta")
    except:
        print("Alignment was empty or malformed. Skipping.")
        pass
    else:
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus = summary_align.dumb_consensus(threshold, "N")
        my_seqs = SeqRecord(
            consensus, id="", description=f"{locus_name}_{threshold*100}%_consensus"
        )
        create_parent_directory(output_file)
        SeqIO.write(my_seqs, output_file, "fasta")



def get_all_tips(fasta: str, tip_list: str = []):
    for seq_record in SeqIO.parse(fasta, "fasta"):
        locus_name = ".".join(fasta.split("/")[-1].split(".")[0:-1])
        name = seq_record.description.replace(f"{locus_name}", "")
        if not name in tip_list:
            tip_list.append(name)
    return tip_list
