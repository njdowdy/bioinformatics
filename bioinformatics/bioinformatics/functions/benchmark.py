import csv
import os
import re
from turtle import st
from typing import Optional
import numpy as np

from bioinformatics.models.phylogeny import PhyloSoftware
from bioinformatics.functions.phylogeny import infer_phylogeny, generate_random_tree
from bioinformatics.functions.file_utils import (
    get_platform_info,
    create_parent_directory,
)
from pyparsing import Opt


def is_partition_file(filename: str, ext_list: list[str] = ["part"]) -> Optional[str]:
    ext_list = "|".join(["." + x + "$" for x in ext_list])
    match = re.search(rf"{ext_list}", filename)
    if match:
        return match.group()


def benchmark_state(file: str, all_files: list[str]) -> str:
    if file + ".bionj" in all_files:
        return "completed"
    elif file + ".ckp.gz" in all_files:
        return "checkpointed"
    elif file + ".log" in all_files:
        return "started"
    else:
        return "queued"


def interrupted_benchmark_cleanup(file: str) -> None:
    os.remove(file + ".log")
    os.remove(file + ".model.gz")


def status_tally(
    alignment_type: Optional[str] = "phy",
    alignment_sets: Optional[str] = "bioinformatics/input/benchmark/sim_data",
) -> None:
    completed = 0
    queued = 0
    started = 0
    other = 0
    total = 0
    todo = []
    for root, dirs, files in os.walk(alignment_sets):
        for filename in [x for x in files if re.search(rf"\.{alignment_type}$", x)]:
            state = benchmark_state(filename, files)
            if state == "completed":
                completed += 1
            elif state == "queued":
                queued += 1
                todo.append(filename)
            elif state == "started":
                started += 1
                todo.append(filename)
            else:
                other += 1
            total += 1
    return {
        "completed": completed,
        "queued": queued,
        "started": started,
        "total": total,
        "remaining": total - completed,
        "todo": todo,
    }


def benchmark_phylogeny(
    builder: PhyloSoftware,
    alignment_type: Optional[str] = "phy",
    alignment_sets: Optional[str] = "bioinformatics/input/benchmark/sim_data",
    use_partition: Optional[bool] = False,
    stdout: Optional[bool] = False,
) -> None:
    platform_info = get_platform_info()
    for root, dirs, files in os.walk(alignment_sets):
        for filename in [
            x
            for x in files
            if re.search(rf"\.{alignment_type}$", x) and not is_partition_file(x)
        ]:
            state = benchmark_state(filename, files)
            if state in ["queued", "started", "checkpointed"]:
                alignment_path = root + "/" + filename
                if state in ["started"]:
                    interrupted_benchmark_cleanup(alignment_path)
                partition_file = None
                if use_partition:
                    part_file = re.sub(
                        r"(\..*?)$", is_partition_file(filename), filename
                    )
                    if part_file in files:
                        partition_file = part_file
                # print(alignment_path)
                infer_phylogeny(alignment_path, builder, partition_file, stdout)


def input_check(range_data: list[float, float, float]) -> bool:
    assert (
        range_data[0] != 0 and range_data[1] != 0
    )  # start and end point should never be zero
    assert (
        range_data[0] <= range_data[1]
    )  # start point should be less than equal to end point
    if range_data[1] - range_data[0] == 0:
        assert (
            range_data[2] == range_data[1]
        )  # if start and end same, step size should be same value
    else:
        assert (range_data[1] - range_data[1]) % range_data[
            2
        ] == 0  # distance to travel evenly divisible by step size
        assert (
            range_data[2] <= range_data[1] - range_data[0]
        )  # step size should be less than distance to travel
    return True


def generate_benchmark_data(
    ntaxa_range: list[int, int, int],
    nsites_range: list[int, int, int],
    processes: Optional[list[str]] = ["bd"],  # TODO: make a enum class to limit options
    output_folder: Optional[str] = "bioinformatics/input/benchmark/sim_data",
    birth_rate_range: Optional[list[float, float, float]] = [0.1, 0.1, 0.0],
    death_rate_range: Optional[list[float, float, float]] = [0.05, 0.05, 0.0],
    stdout: Optional[bool] = False,
) -> None:
    # fix input ranges to get proper end point (zero indexing problem)
    ntaxa_range, nsites_range, birth_rate_range, death_rate_range = list(
        map(
            lambda x: [x[0], x[1] + x[2], x[2]],
            [ntaxa_range, nsites_range, birth_rate_range, death_rate_range],
        )
    )
    # change a step size of 0 to x[1] for np.arange
    ntaxa_range, nsites_range, birth_rate_range, death_rate_range = list(
        map(
            lambda x: [x[0], x[1] + x[1], x[1]] if x[2] == 0 else [x[0], x[1], x[2]],
            [ntaxa_range, nsites_range, birth_rate_range, death_rate_range],
        )
    )
    verify_inputs = map(
        input_check,
        [ntaxa_range, nsites_range, birth_rate_range, death_rate_range],
    )
    if all(list(verify_inputs)):  # verify passed data
        for process in processes:
            for ntaxa in [
                x for x in np.arange(ntaxa_range[0], ntaxa_range[1], ntaxa_range[2])
            ]:
                for nsites in [
                    y
                    for y in np.arange(
                        nsites_range[0], nsites_range[1], nsites_range[2]
                    )
                ]:
                    for birth_rate in [
                        z
                        for z in np.arange(
                            birth_rate_range[0],
                            birth_rate_range[1],
                            birth_rate_range[2],
                        )
                    ]:
                        for death_rate in [
                            t
                            for t in np.arange(
                                death_rate_range[0],
                                death_rate_range[1],
                                death_rate_range[2],
                            )
                        ]:
                            generate_random_tree(
                                ntaxa=ntaxa,
                                nsites=nsites,
                                process=process,
                                output_folder=output_folder,
                                birth_rate=birth_rate,
                                death_rate=death_rate,
                                stdout=stdout,
                            )
    else:
        raise Exception(ValueError)


def extract_text_matches(result: dict[str, re.Match]) -> dict[str, str]:
    for key, val in dict.items(result):
        if val:
            val = val.group(1)
            result[f"{key}"] = val
        else:
            result[f"{key}"] = None
    return result


def parse_log_files(
    alignment_type: Optional[str] = "phy",
    alignment_sets: Optional[str] = "bioinformatics/input/benchmark/sim_data",
    results_file: Optional[
        str
    ] = "bioinformatics/output/benchmark/benchmark_results.csv",
) -> None:
    fields = [
        "ntaxa",
        "nsites",
        "process",
        "birth_rate",
        "death_rate",
        "fast_ml_tree_time",
        "model_cpu_time",
        "model_wall_time",
        "tree_cpu_time",
        "tree_wall_time",
        "model_ram_used",
        "tree_ram_used",
        "host_info",
        "threads",
        "best_fit_model",
    ]
    results = []
    for root, dirs, files in os.walk(alignment_sets):
        for logfile in [x for x in files if re.search(rf"\.{alignment_type}\.log$", x)]:
            result = {}
            result["ntaxa"] = re.search(r"_([0-9]*?)_taxa", logfile)
            result["nsites"] = re.search(r"_([0-9]*?)_sites", logfile)
            result["process"] = re.search(r"benchmark_data_(.*?)_", logfile)
            result["birth_rate"] = re.search(r"br_(.*?)_dr", logfile)
            result["death_rate"] = re.search(r"dr_(.*?)_", logfile)
            with open(root + "/" + logfile, "r") as log:
                text = "".join(log.readlines())
                result["host_info"] = re.search(r"Host:\s+(.*?)\n", text)
                result["threads"] = re.search(r" - ([0-9]*?) threads", text)
                result["fast_ml_tree_time"] = re.search(
                    r"fast ML tree search: (.*?) sec", text
                )
                result["best_fit_model"] = re.search(
                    r"Best-fit model: (.*?) chosen", text
                )
                result["model_cpu_time"] = re.search(
                    r"CPU time for ModelFinder: (.*?) sec", text
                )
                result["model_wall_time"] = re.search(
                    r"Wall-clock time for ModelFinder: (.*?) sec", text
                )
                result["tree_cpu_time"] = re.search(
                    r"CPU time used for tree search: (.*?) sec", text
                )
                result["tree_wall_time"] = re.search(
                    r"Wall-clock time used for tree search: (.*?) sec", text
                )
                result["model_ram_used"] = re.search(
                    r"NOTE: ModelFinder requires ([0-9]*?) MB RAM!", text
                )
                result["tree_ram_used"] = re.search(r"NOTE: ([0-9]*?) MB RAM", text)
            results.append(extract_text_matches(result))
    create_parent_directory(results_file)
    with open(results_file, "w") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
