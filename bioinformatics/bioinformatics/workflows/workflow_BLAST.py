from bioinformatics.functions.tasks import (
    subset_fasta_alignments,
    generate_consensus_seqs,
    pad_alignments,
)
from bioinformatics.functions.file_utils import combine_files
from bioinformatics.functions.alignment import perform_alignment
from bioinformatics.functions.blast import create_blast_db

# step 1
subset_fasta_alignments()
# step 2
pad_alignments()
# step 3
perform_alignment()
# step 4 (optional) - trim alignments

# step 5
generate_consensus_seqs()
# step 6
combine_files()
# step 7 - generate BLAST databases
create_blast_db()
