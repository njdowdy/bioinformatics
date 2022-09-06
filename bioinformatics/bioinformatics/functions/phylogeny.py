from typing import Optional
from bioinformatics.models.phylogeny import PhyloSoftware


from bioinformatics.functions.file_utils import (
    generate_process,
    create_parent_directory,
)


def infer_phylogeny(
    alignment: str,
    builder: PhyloSoftware,
    partition_file: Optional[str] = None,
    stdout: bool = False,
) -> None:
    # alignment = 'bioinformatics/input/test_files/example.phy'
    # builder = PhyloSoftware("iqtree2")
    # partition_file = 'bioinformatics/input/test_files/example.nex'
    cmd = []
    if builder.name == "iqtree2":
        src = "bioinformatics/src/phylogenetics/iqtree2_v2.2.0"
        in_param = "-s"
        cmd.extend([src, in_param, alignment])
        if partition_file:
            partition_param = "-p"
            cmd.extend([partition_param, partition_file])
        cmd.extend(["-T", "6"])
        cmd.extend(["-m", "MFP"])
        # cmd.extend(["-mset", "JC,F81,HKY,TIM,GTR"])
    elif builder.name == "raxmlng":
        src = "bioinformatics/src/phylogenetics/raxml-ng_v1.1.0"
        in_param = "--msa"
        cmd.extend([src, in_param, alignment])
        if partition_file:
            partition_param = "--model"
            cmd.extend([partition_param, partition_file])
        else:
            partition_param = "--model"
            partition_file = "GTR+G"
            cmd.extend([partition_param, partition_file])
    else:
        raise Exception("Phylogeny software choice not understood.")
    print(cmd)
    generate_process(cmd, stdout=True)


def generate_random_tree(
    ntaxa: int,
    nsites: int,
    process: Optional[str] = "bd",
    output_folder: Optional[str] = "bioinformatics/input/benchmark/sim_data/",
    birth_rate: Optional[float] = 0.1,
    death_rate: Optional[float] = 0.05,
    stdout: Optional[bool] = False,
):
    # more info and options: http://www.iqtree.org/doc/AliSim
    cmd = []
    src = "bioinformatics/src/phylogenetics/iqtree2_v2.2.0"
    cmd.extend(
        [
            src,
            "--alisim",
        ]
    )
    if process == "bd":
        output = f"{output_folder}/benchmark_data_{process}_br_{birth_rate}_dr_{death_rate}_{ntaxa}_taxa_{nsites}_sites"
        create_parent_directory(output)
        cmd.extend(
            [
                f"{output}",
                "-t",
                f"RANDOM{{{process}{{{birth_rate}/{death_rate}}}/{ntaxa}}}",
            ]
        )
    elif process in ["yh", "u", "cat", "bal"]:
        output = f"{output_folder}/benchmark_data_{process}_{ntaxa}_taxa_{nsites}_sites"
        create_parent_directory(output)
        cmd.extend(
            [
                f"{output}",
                "-t",
                f"RANDOM{{{process}/{ntaxa}}}",
            ]
        )
    else:
        raise Exception(ValueError)
    cmd.extend(["--length", f"{nsites}"])

    generate_process(cmd, stdout=stdout)
