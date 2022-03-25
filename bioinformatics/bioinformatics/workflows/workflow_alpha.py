import os
from typing import List, Optional
from bioinformatics.functions.ingest import filter_fasta_by_label


def subset_fasta_alignments(
    input_folder: str, label_filter: List[str], backup_filter: Optional[List[str]]
):
    for root, dirs, files in os.walk(input_folder):
        for filename in files:
            print(os.path.join(root, filename))
            # filter_fasta_by_label(
            #     fasta=os.path.join(root, filename),
            #     label_filter=label_filter,
            #     backup_filter=backup_filter,
            # )
