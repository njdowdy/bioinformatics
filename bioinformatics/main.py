import os
from bioinformatics.workflows.workflow_alpha import (
    subset_fasta_alignments,
    generate_consensus_seqs,
    get_all_tips_labels,
    pad_alignments,
)
from bioinformatics.utils.zfmk1 import filters as zfmk1
from bioinformatics.utils.lep1 import filters as lep1

input_folder = "bioinformatics/input/fasta/LEP1"
folder2 = "bioinformatics/output/taxon_filtered_alignments/LEP1"

primary_filter = lep1.primary_filter
secondary_filter = lep1.secondary_filter


def main():
    # get_all_tips_labels(
    #     input_folder,
    #     f"bioinformatics/utils/{input_folder.split('/')[-1]}/all_labels.py",
    # )
    # subset_fasta_alignments(input_folder, primary_filter, secondary_filter)
    # TODO: check that they are aligned
    # generate_consensus_seqs(folder2, 0.5)
    pad_alignments(input_folder)


if __name__ == "__main__":
    main()
