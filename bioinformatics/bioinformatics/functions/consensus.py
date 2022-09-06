import re
from typing import Optional, List
from Bio import SeqIO
from Bio.SeqIO import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo
from bioinformatics.functions.file_utils import (
    parse_output_file_name,
    create_parent_directory,
)


def generate_consensus(fasta: str, threshold: Optional[float] = 0.5):
    # fasta = 'bioinformatics/output/taxon_filtered_alignments/ZFMK1/EOG090R0A51_3.final.out.fas'
    output_file = parse_output_file_name(fasta, 'output/consensus_sequences').replace(".out.out", ".out")
    # subfolders_to_replace = "/".join(
    #     output_file.split("/")[output_file.split("/").index("output") + 1 : -2]
    # )
    # output_file = output_file.replace(f"{subfolders_to_replace}", "consensus_sequences")
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

