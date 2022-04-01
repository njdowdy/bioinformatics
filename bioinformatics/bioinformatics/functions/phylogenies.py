import subprocess
from typing import Optional
from bioinformatics.models.phylogeny import PhyloSoftware


def infer_phylogeny(
    alignment: str, builder: PhyloSoftware, partition_file: Optional[str] = None
) -> None:
    # alignment = 'bioinformatics/input/test_files/example.phy'
    # builder = PhyloSoftware("iqtree2")
    # partition_file = 'bioinformatics/input/test_files/example.nex'
    cmd = []
    if builder.name == "iqtree2":
        src = "bioinformatics/src/phylogenetics/iqtree2_v2.1.2"
        in_param = "-s"
        cmd.extend([src, in_param, alignment])
        if partition_file:
            partition_param = "-p"
            cmd.extend([partition_param, partition_file])
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
    subprocess.run(
        cmd,
    )
