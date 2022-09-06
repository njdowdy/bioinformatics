from distutils.command.build import build
from bioinformatics.functions.tasks import (
    generate_consensus_seqs,
    create_blast_dbs,
    replace_char_in_filenames,
    apply_single_fasta_blast,
    apply_fasta_clean_newlines,
)
from bioinformatics.functions.tasks import (
    pad_alignments,
    perform_alignments,
    fix_dna_alignments,
    apply_replace_ambiguous_chars,
)
from bioinformatics.functions.file_utils import copy_directory
from bioinformatics.models.alignment import (
    AlignmentOutputFormat,
    AlignmentSoftware,
)
from bioinformatics.models.phylogeny import PhyloSoftware
from bioinformatics.functions.partition import amas_create_supermatrix
from bioinformatics.functions.phylo_tips import drop_taxa
from bioinformatics.functions.phylogeny import infer_phylogeny

input_folder = "bioinformatics/input/fasta/NOC1_PHASE2/"
reference_seqs = "bioinformatics/input/reference_seqs/"

# clean reference filenames
replace_char_in_filenames(reference_seqs, " ")

# generate consensus sequences
generate_consensus_seqs(input_folder)

# generate BLAST database for each reference
create_blast_dbs(reference_seqs)  # ~683 loci mapped; ~133 unmapped

# step 8: BLAST each locus > parse top hit as extracted region
extraction_folder = "bioinformatics/output/seq_extractions/NOC1_PHASE2/"
copy_directory(input_folder, extraction_folder)
input_folder = "bioinformatics/output/consensus_sequences/NOC1_PHASE2/"
ref_blast_db = "bioinformatics/output/blastdb/genome_spilosoma_lubricepidum/"
apply_single_fasta_blast(
    input_folder=input_folder,
    ref_blast_db=ref_blast_db,
    extraction_folder=extraction_folder,
    dcmegablast=False,
)
apply_fasta_clean_newlines(extraction_folder)

# step 10: Align FASTA file
input_folder = "bioinformatics/output/seq_extractions/NOC1_PHASE2/"
padded_input_folder = input_folder + "padded/"

pad_alignments(input_folder)
perform_alignments(
    padded_input_folder, AlignmentSoftware("muscle3"), AlignmentOutputFormat("fasta")
)
orf_input_folder = "bioinformatics/output/alignments/muscle3/NOC1_PHASE2/padded/"
fix_dna_alignments(orf_input_folder)

# step 11: concatenate supermatrix
input_folder = "bioinformatics/output/alignments/muscle3/NOC1_PHASE2/ORF1"
supermatrix_outfile = "bioinformatics/output/supermatrix/NOC1_PHASE2/alignment.phy"
partition_outfile = "bioinformatics/output/supermatrix/NOC1_PHASE2/alignment.part"
amas_create_supermatrix(
    input_folder=input_folder,
    input_format="fasta",
    supermatrix_outfile=supermatrix_outfile,
    partition_outfile=partition_outfile,
    cores=4,
    supermatrix_format="phylip",
    partition_format="nexus",
    codons="123",
)
apply_replace_ambiguous_chars(supermatrix_outfile, ["?"], "N")

# step 12: prune taxa (~30 taxa)
input_file = "bioinformatics/output/supermatrix/NOC1_PHASE2/alignment.phy"
output_file = "bioinformatics/output/supermatrix/NOC1_PHASE2/alignment_sub35.phy"
taxa_to_retain = [
    "_Spilosoma_lubricepidum_BLAST_genome_extracted",
    "Spilosoma_vagans",
    "Spilarctia_nydia_tienmushanica",
    "Hyphantria_cunea",
    "Haploa_confusa",
    "Virbia_nigricans",
    "Hypercompe_laeta",
    "Acyphas_chionitis",
    "Lymantria_dispar",
]
drop_taxa(
    input_file=input_file,
    output_file=output_file,
    final_taxa_count=30,
    taxa_to_retain=taxa_to_retain,
    taxa_to_remove=[],
)

# step 13: basic IQTREE2 run
input_file = "bioinformatics/output/supermatrix/NOC1_PHASE2/alignment_sub35.phy"
infer_phylogeny(
    alignment=input_file,
    partition_file=partition_outfile,
    builder=PhyloSoftware("iqtree2"),
)
