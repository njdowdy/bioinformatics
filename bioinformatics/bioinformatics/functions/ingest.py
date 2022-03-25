import os
import re
from typing import Optional, List
from Bio import SeqIO

def filter_fasta_by_label(
    fasta: str, label_filter: List[str], backup_filter: Optional[List[str]]
):
    suffix = re.search(r"\..*?$", fasta).group()
    output_file = re.sub(rf"{suffix}$", f"_out{suffix}", fasta)
    output_file = re.sub(r"/input/", "/output/", output_file)
    os.makedirs("/".join(output_file.split("/")[0:-1]), exist_ok=True)
    with open(output_file, "w") as f:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            lab_found = False
            for substring in label_filter:
                if substring in seq_record.description:
                    f.write(f">{seq_record.id}\n{seq_record.seq}\n")
                    lab_found = True
            if not lab_found:
                for substring in backup_filter:
                    if substring in seq_record.description:
                        f.write(f">{seq_record.id}\n{seq_record.seq}\n")
    pass
