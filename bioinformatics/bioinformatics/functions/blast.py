import subprocess
from typing import List, Optional
from bioinformatics.functions.file_utils import (
    replace_parent_directory,
    create_parent_directory,
)


def create_blast_db(
    consensus_seq_set: str, blast_db_out: Optional[str], data_type: str = "nucl"
) -> None:
    if not blast_db_out:
        blast_db_out = replace_parent_directory(
            consensus_seq_set, "output/consensus_sequences", "output/blastdb"
        )
    create_parent_directory(blast_db_out)
    cmd = [
        "bioinformatics/src/blast/ncbi-blast-2.13.0+/bin/makeblastdb",
        "-in",
        consensus_seq_set,
        "-parse_seqids",
        "-dbtype",
        data_type,
        "-out",
        blast_db_out,
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
    if not translated:
        blast_program = "bioinformatics/src/blast/ncbi-blast-2.13.0+/bin/blastn"
    else:
        blast_program = "bioinformatics/src/blast/ncbi-blast-2.13.0+/bin/tblastn"
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
        cmd.extend(["-task", "dc-megablast"])
    subprocess.run(
        cmd,
    )


def pair_all_blast_dbs(blast_dbs: List[str]) -> List[tuple[str]]:
    result = [(a, b) for idx, a in enumerate(blast_dbs) for b in blast_dbs[idx + 1 :]]
    return result
