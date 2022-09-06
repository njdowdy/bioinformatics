from bioinformatics.functions.tasks import (
    pad_alignments,
    perform_alignments,
    trim_alignments,
    fix_dna_alignments,
    apply_replace_ambiguous_chars,
)
from bioinformatics.models.alignment import (
    AlignmentOutputFormat,
    AlignmentSoftware,
)
from bioinformatics.models.trim import (
    TrimSoftware,
)
from bioinformatics.functions.file_utils import list_files

input_folder = "bioinformatics/input/fasta/LEP1/"
padded_input_folder = input_folder + "padded/"
orf_input_folder = "bioinformatics/output/alignments/muscle3/LEP1/padded/"
trimming_input_folder = "bioinformatics/output/alignments/muscle3/LEP1/ORF1/"
trimmed_output_folder = "bioinformatics/output/alignments/muscle3/LEP1/trimmed/bmge/"

# step 0
list_files(input_folder)
# step 1
pad_alignments(input_folder)
# step 2
perform_alignments(
    padded_input_folder, AlignmentSoftware("muscle3"), AlignmentOutputFormat("fasta")
)
# step 3
fix_dna_alignments(orf_input_folder)
# step 4 (optional) - trim alignments
trim_alignments(trimming_input_folder, TrimSoftware("trimal"))
trim_alignments(trimming_input_folder, TrimSoftware("trimal"), method="gappyout")
trim_alignments(trimming_input_folder, TrimSoftware("trimal"), method="strictplus")
trim_alignments(trimming_input_folder, TrimSoftware("trimal"), method="strict")
trim_alignments(
    trimming_input_folder, TrimSoftware("bmge")
)  # this is preferred because CODON is available, which preserves codon ORF
trim_alignments(trimming_input_folder, TrimSoftware("clipkit"))
# step 5 (optional) - clean up ambiguous chars
apply_replace_ambiguous_chars(trimmed_output_folder, ["X"], "?")
# step 6 (optional) - code gap chars
