import os
import re
from tqdm import tqdm
from typing import List, Optional
from bioinformatics.functions.ingest import (
    filter_fasta_by_label,
    get_all_tips,
)
from bioinformatics.functions.blast import (
    create_blast_db,
    single_fasta_blast,
    inject_blast_result,
)
from bioinformatics.functions.consensus import generate_consensus
from bioinformatics.functions.align import (
    pad_alignment,
    perform_alignment,
)
from bioinformatics.functions.ingest import get_molecular_data_type
from bioinformatics.functions.clean import (
    replace_char_in_filename,
    fasta_clean_newlines,
    replace_ambiguous_chars,
)
from bioinformatics.models.alignment import (
    AlignmentOutputFormat,
    AlignmentSoftware,
)
from bioinformatics.models.trim import TrimSoftware
from bioinformatics.functions.trim import trim_alignment
from bioinformatics.functions.orf import fix_dna_alignment
from bioinformatics.functions.file_utils import suffix_parser


def subset_fasta_alignments(
    input_folder: str, primary_filter: List[str], secondary_filter: Optional[List[str]]
):
    for root, dirs, files in os.walk(input_folder):
        for filename in tqdm(files):
            filter_fasta_by_label(
                fasta=os.path.join(root, filename),
                primary_filter=primary_filter,
                secondary_filter=secondary_filter,
            )


def generate_consensus_seqs(input_folder: str, threshold: Optional[float] = 0.5):
    for root, dirs, files in os.walk(input_folder):
        for filename in tqdm(files):
            generate_consensus(fasta=os.path.join(root, filename), threshold=threshold)


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
        for filename in tqdm(files):
            pad_alignment(os.path.join(root, filename))


def perform_alignments(
    input_folder: str,
    aligner: AlignmentSoftware,
    output_format: AlignmentOutputFormat,
    iterations: Optional[int] = None,
    stdout: bool = False,
) -> None:
    for root, dirs, files in os.walk(input_folder):
        for filename in tqdm(files):
            perform_alignment(
                os.path.join(root, filename),
                aligner,
                output_format,
                iterations,
                stdout,
            )


def trim_alignments(
    input_folder: str,
    trimmer: TrimSoftware,
    method: Optional[str] = None,
    stdout: bool = False,
) -> None:
    for root, dirs, files in os.walk(input_folder):
        for filename in tqdm(files):
            trim_alignment(os.path.join(root, filename), trimmer, method, stdout)


def fix_dna_alignments(input_folder: str, best_orf: Optional[int] = None) -> None:
    for root, dirs, files in os.walk(input_folder):
        for filename in tqdm(files):
            fix_dna_alignment(os.path.join(root, filename), best_orf)


def apply_replace_ambiguous_chars(
    input_folder: str, search_chars: list[str], replace_char: Optional[str] = "?"
):
    for root, dirs, files in os.walk(input_folder):
        for filename in tqdm(files):
            replace_ambiguous_chars(
                os.path.join(root, filename), search_chars, replace_char
            )


def replace_char_in_filenames(
    input_folder: str, search_char: str, replace_char: Optional[str] = "_"
) -> None:
    for root, dirs, files in os.walk(input_folder):
        for filename in tqdm(files):
            replace_char_in_filename(
                os.path.join(root, filename),
                search_char=" ",
                replacement_char=replace_char,
            )


def create_blast_dbs(
    input_folder: str,
    blast_db_out: Optional[str] = None,
) -> None:
    for root, dirs, files in os.walk(input_folder):
        for filename in tqdm(files):
            try:
                data_type = get_molecular_data_type(os.path.join(root, filename))
                create_blast_db(os.path.join(root, filename), blast_db_out, data_type)
            except:
                pass


def apply_single_fasta_blast(
    input_folder: str,
    ref_blast_db: str,
    extraction_folder: str,
    translated: Optional[bool] = False,
    dcmegablast: Optional[bool] = False,
) -> None:
    ref_blast_db_tag = re.sub(r"/$", "", ref_blast_db).split("/")[-1]
    ref_blast_db = ref_blast_db + f"{ref_blast_db_tag}"
    for root, dirs, files in os.walk(input_folder):
        for filename in tqdm(files):
            if dcmegablast:
                blast_result_prefix = "bioinformatics/output/blast_hits/dc-megablast/"
            else:
                blast_result_prefix = "bioinformatics/output/blast_hits/blastn/"
            query_set_tag = re.sub(r"/$", "", root).split("/")[-1]
            query_tag = re.sub(rf"\.{suffix_parser(filename)}", "", filename)
            blast_result_out = f"{blast_result_prefix}{query_set_tag}_VS_{ref_blast_db_tag}/{query_tag}"
            single_fasta_blast(
                query=os.path.join(root, filename),
                ref_db=ref_blast_db,
                blast_result_out=blast_result_out,
                translated=translated,
                dcmegablast=dcmegablast,
            )
            inject_blast_result(
                input_file=os.path.join(extraction_folder, filename),
                blast_hit_file=blast_result_out,
            )


def apply_fasta_clean_newlines(
    input_folder: str,
) -> None:
    for root, dirs, files in os.walk(input_folder):
        for filename in tqdm(files):
            fasta_clean_newlines(os.path.join(root, filename))
