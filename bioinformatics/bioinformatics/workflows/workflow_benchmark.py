from bioinformatics.functions.benchmark import (
    generate_benchmark_data,
    benchmark_phylogeny,
    parse_log_files,
)
from bioinformatics.models.phylogeny import PhyloSoftware

ntaxa_range = [10, 200, 40]
nsites_range = [10000, 460000, 50000]
processes = ["bd"]
birth_rate_range = [0.1, 0.1, 0]
death_rate_range = [0.05, 0.05, 0]
generate_benchmark_data(
    ntaxa_range=ntaxa_range,
    nsites_range=nsites_range,
    processes=processes,
    birth_rate_range=birth_rate_range,
    death_rate_range=death_rate_range,
    stdout=False,
)
benchmark_phylogeny(builder=PhyloSoftware("iqtree2"), use_partition=False, stdout=True)
parse_log_files()
