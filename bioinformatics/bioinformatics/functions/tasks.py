import os
from typing import List, Optional
from bioinformatics.functions.ingest import (
    filter_fasta_by_label,
    generate_consensus,
    get_all_tips,
)
from bioinformatics.functions.alignment import pad_alignment


def subset_fasta_alignments(
    input_folder: str, primary_filter: List[str], secondary_filter: Optional[List[str]]
):
    for root, dirs, files in os.walk(input_folder):
        for filename in files:
            filter_fasta_by_label(
                fasta=os.path.join(root, filename),
                primary_filter=primary_filter,
                secondary_filter=secondary_filter,
            )


def generate_consensus_seqs(input_folder: str, threshold: Optional[float]):
    for root, dirs, files in os.walk(input_folder):
        for filename in files:
            generate_consensus(fasta=os.path.join(root, filename))


def get_all_tips_labels(input_folder: str, output_file: str):
    tip_list_complete = []
    for root, dirs, files in os.walk(input_folder):
        for filename in files:
            tip_list = get_all_tips(
                fasta=os.path.join(root, filename), tip_list=tip_list_complete
            )
            tip_list_complete = list(set(tip_list + tip_list_complete))
    os.makedirs("/".join(output_file.split("/")[0:-1]), exist_ok=True)
    with (open(output_file, "w")) as f:
        f.write(f"{set(tip_list_complete)}")


def pad_alignments(input_folder: str) -> None:
    for root, dirs, files in os.walk(input_folder):
        for filename in files:
            pad_alignment(os.path.join(root, filename))
