import os
from bioinformatics.workflows.workflow_alpha import subset_fasta_alignments

input_folder = "bioinformatics_tools/input"
label_filter = ["Erebidae_"]
backup_filter = [
    "Noctuidae_",
    "Notodontidae_",
    "Nolidae_",
    "Euteliidae_",
]


def main():
    subset_fasta_alignments(input_folder, label_filter, backup_filter)


if __name__ == "main":
    main()
