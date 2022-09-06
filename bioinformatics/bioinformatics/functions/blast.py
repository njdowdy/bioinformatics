import re
import subprocess
from typing import List, Optional
from bioinformatics.functions.file_utils import (
    replace_parent_directory,
    create_parent_directory,
)


def get_blast_program(translated: bool = False) -> str:
    if not translated:
        blast_program = "bioinformatics/src/blast/ncbi-blast-2.13.0+/bin/blastn"
    else:
        blast_program = "bioinformatics/src/blast/ncbi-blast-2.13.0+/bin/tblastn"
    return blast_program


def apply_dcmegablast(cmd: subprocess) -> subprocess:
    cmd.extend(["-task", "dc-megablast"])
    return cmd


def create_blast_db(
    consensus_seq_set: str,
    blast_db_out: Optional[str] = None,
    data_type: Optional[str] = "nucl",
) -> None:
    if not blast_db_out:
        consensus_location = "/".join(consensus_seq_set.split("/")[1:-1])
        blast_db_out = "".join(
            replace_parent_directory(
                consensus_seq_set, consensus_location, "output/blastdb"
            ).split(".")[0:-1]
        ).replace(" ", "_")
    create_parent_directory(blast_db_out + "/")
    blast_db_out = blast_db_out + f"/{blast_db_out.split('/')[-1]}"
    cmd = [
        "bioinformatics/src/blast/ncbi-blast-2.13.0+/bin/makeblastdb",
        "-in",
        consensus_seq_set,
        "-parse_seqids",
        "-dbtype",
        data_type,
        "-out",
        blast_db_out,
        "-title",
        blast_db_out.split("/")[-1],
    ]
    subprocess.run(
        cmd,
    )


def reciprocal_blastn(
    ref_db: str,
    query_db: str,
    blast_result_out: str,
    translated: bool = False,
    megablast: bool = False,
) -> None:
    # ../ncbi-blast-2.12.0/bin/blastn -db ./Phase1_BLAST -query ./Phase1_Consensus_Sequences_MajorityRule.fasta -out ../output/blast_results/p1_selfcheck -outfmt 5
    if ref_db == query_db:
        blast_result_out + "_selfcheck"
    create_parent_directory(blast_result_out)
    blast_program = get_blast_program(translated)
    cmd = [
        blast_program,
        "-db",
        ref_db,
        "-query",
        query_db,
        "-out",
        blast_result_out,
        "-outfmt",
        "5",
    ]
    if megablast:
        cmd = apply_dcmegablast(cmd)
    subprocess.run(
        cmd,
    )


def pair_all_blast_dbs(blast_dbs: List[str]) -> List[tuple[str]]:
    result = [(a, b) for idx, a in enumerate(blast_dbs) for b in blast_dbs[idx + 1 :]]
    return result


def single_fasta_blast(
    query: str,
    ref_db: str,
    blast_result_out: str,
    translated: Optional[bool] = False,
    dcmegablast: Optional[bool] = False,
):
    blast_program = get_blast_program(translated)
    create_parent_directory(blast_result_out)
    cmd = [
        blast_program,
        "-db",
        ref_db,
        "-query",
        query,
        "-out",
        blast_result_out,
        "-outfmt",
        "5",
    ]
    if dcmegablast:
        cmd = apply_dcmegablast(cmd)
    subprocess.run(
        cmd,
    )


def parse_blast_result(blast_hit_file: str) -> tuple[str, str]:
    with open(blast_hit_file) as f:
        hit_data = f.read()
    seq = ""
    lab = ""
    first_hit_scoring_pair = re.search(r"<Hsp_hseq>(.*?)</Hsp_hseq>", hit_data)
    if first_hit_scoring_pair:
        seq = first_hit_scoring_pair.group(1)
    label_data = re.search(r"<BlastOutput_db>(.*?)</BlastOutput_db>", hit_data)
    if label_data:
        lab = label_data.group(1).split("/")[-1].replace("genome_", "").capitalize()
        lab = lab + "_BLAST_genome_extracted"
    return (lab, seq)


def inject_blast_result(
    input_file: str, blast_hit_file: str, output_file: Optional[str] = None
) -> None:
    if not output_file:
        output_file = input_file.replace(".out", "")
    create_parent_directory(output_file)
    lab, seq = parse_blast_result(blast_hit_file)
    if lab and seq:
        with open(output_file, "a") as f:
            f.write(f">{lab}\n{seq}\n")
