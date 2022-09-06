import csv
import random
import re
from os import path
import pandas as pd
from typing import Optional
from bioinformatics.functions.file_utils import (
    inject_prefix_suffix,
    create_parent_directory,
)


def form_regex(
    replace_chars: list[str], prefix: Optional[str] = "", suffix: Optional[str] = ""
) -> str:
    regex = "|".join(
        [
            f"{prefix}{y}{suffix}"
            for y in [
                x.replace("(", "\(").replace(")", "\)".replace(".", "\."))
                for x in replace_chars
            ]
        ]
    )
    return regex


def read_new_tips(
    map_file: str,
    search_col_name: str,
    replace_col_name: str,
    replace_chars: Optional[list[str]] = ["(", ")"],
) -> dict[str, str]:
    regex = form_regex(replace_chars)
    df = pd.read_csv(map_file)
    if not df[f"{search_col_name}"].is_unique:
        raise ValueError("Search Column Contains Non-Unique Labels!")
    name_maps = dict(
        zip(
            # df[f"{search_col_name}"].str.replace(r"\.", "_", regex=True),
            df[f"{search_col_name}"],
            df[f"{replace_col_name}"].str.replace(rf"{regex}", "_", regex=True),
        )
    )
    return name_maps


def replace_tips(
    in_file: str,
    map_file: str,
    search_col_name: Optional[str] = "original_tip",
    replace_col_name: Optional[str] = "new_tip",
) -> None:
    # in_file = 'bioinformatics/input/phylo_tips/tree.tre'
    # map_file = 'bioinformatics/input/phylo_tips/data3.csv'
    # search_col_name = 'Taxon.Code'
    # replace_col_name = 'Species.nameNEW'
    # replace_tips(in_file, map_file, search_col_name, replace_col_name)
    name_maps = read_new_tips(map_file, f"{search_col_name}", f"{replace_col_name}")
    tree = open(in_file, "r").read()
    # tree = re.sub(
    #     rf"({'|'.join(name_maps.keys())})", lambda m: rf"{name_maps[m.group(1)]}", tree
    # )
    for k, v in name_maps.items():
        tree = tree.replace(k, v)
    out_file = inject_prefix_suffix(in_file, "", "_new_tips")
    create_parent_directory(out_file)
    with open(out_file, "w") as f:
        f.write(tree)


def get_taxa(in_file: str) -> list[str]:
    taxa = []
    with open(in_file) as f:
        lines = f.readlines()[1:]
        for line in lines:
            taxon = re.search(r"^[^ ]* ", line)
            if taxon:
                taxa.append(taxon.group(0))
    return taxa


def get_taxa_from_phylo(in_file: str) -> Optional[list[str]]:
    with open(in_file) as f:
        tree = f.readlines()
    if len(tree) > 0:
        taxa = re.findall(r"[\(|,](.*?):", tree[0])
        taxa = [re.sub(r"\(|\)|,", "", x) for x in taxa]
    else:
        taxa = None
    return taxa


def epithet_extractor(
    taxa: list[str],
    regex: Optional[str] = "[A-Z][a-z]+_[a-z]+(?!.*[A-Z][a-z])",
    replace_chars: Optional[list[str]] = ["cf", "nr"],
) -> list[str]:
    replace = form_regex(replace_chars, "_")
    new_taxa = [re.findall(rf"{regex}", x) for x in taxa]
    new_taxa = [
        re.sub(rf"{replace}", "_", x[0]).replace("_", " ") if x else ""
        for x in new_taxa
    ]
    new_taxa = [re.sub(r"sp$", "sp.", x) for x in new_taxa]
    ids = [re.findall(rf"^(.*?)_", x) for x in taxa]
    ids = [f" [{x[0]}]" if x else "" for x in ids]
    new_labels = []
    for count, value in enumerate(new_taxa):
        new_labels.append(value + ids[count])
    return new_labels


def write_tip_map(
    taxa: list[str],
    out_file: str,
    in_file: Optional[str] = None,
    header: Optional[str] = "original_tip",
) -> None:
    if in_file:
        if path.exists(in_file):
            csv_data = pd.read_csv(in_file)
    elif path.exists(out_file):
        csv_data = pd.read_csv(out_file)
    else:
        csv_data = pd.DataFrame()
    csv_data[f"{header}"] = taxa
    csv_data.to_csv(out_file, index=False)


def get_taxa_count(in_file: str) -> int:
    with open(in_file) as f:
        first_line = f.readline().split(" ")
        if first_line:
            return int(first_line[0])


def get_site_count(in_file: str) -> int:
    with open(in_file) as f:
        first_line = f.readline().split(" ")
        if first_line:
            return int(first_line[1])


special_chars = ["[", "]", "(", ")", "{", "}", "*", "+", "?", "|", "^", "$", ".", "\\"]


def drop_taxa(
    input_file: str,
    final_taxa_count: int,
    output_file: Optional[str] = "",
    taxa_to_retain: Optional[list[str]] = [],
    taxa_to_remove: Optional[list[str]] = [],
) -> None:
    # get taxa
    taxa = get_taxa(input_file)
    assert (
        len(taxa) >= final_taxa_count
    ), "Desired final alignment size is larger than current alignment size."
    assert len(taxa_to_retain) <= len(
        taxa
    ), "Cannot retain more taxa than present in alignment."
    assert len(taxa_to_remove) <= len(
        taxa
    ), "Cannot remove more taxa than present in alignment."
    assert (
        len(taxa_to_retain) <= final_taxa_count
    ), "Cannot retain more than the desired final alignment size."
    assert len(taxa) - len(taxa_to_remove) >= final_taxa_count - len(
        taxa_to_retain
    ), "Cannot drop that many taxa and meet desired final alignment size."

    # create drop pool
    drop_pool = len(taxa) - len(taxa_to_remove)
    # create random sample to retain
    random_sample = final_taxa_count - len(taxa_to_retain)
    # verify
    assert (
        random_sample <= drop_pool
    ), "Cannot drop/retain that many taxa and meet desired final alignment size."
    with open(input_file) as f:
        lines = f.readlines()[1:]
        retained_lines = []
        # retain taxa to retain
        if taxa_to_retain:
            taxa_to_retain_clean = [
                "".join([f"\{y}" if y in special_chars else y for y in x])
                for x in taxa_to_retain
            ]
            retain_regex = "|".join(taxa_to_retain_clean)
            _ = [
                retained_lines.append(line)
                for line in lines
                if re.search(rf"{retain_regex}", line)
            ]
        # drop taxa to remove
        if taxa_to_remove:
            taxa_to_remove_clean = [
                "".join([f"\{y}" if y in special_chars else y for y in x])
                for x in taxa_to_remove
            ]
            remove_regex = "|".join(taxa_to_remove_clean)
            if taxa_to_retain:
                remove_regex = remove_regex + "|" + retain_regex
            new_lines = [
                line for line in lines if not re.search(rf"{remove_regex}", line)
            ]
        elif taxa_to_retain:
            new_lines = [
                line for line in lines if not re.search(rf"{retain_regex}", line)
            ]
        else:
            new_lines = lines
        # randomly sample lines from pool
        if random_sample > 0:
            lines_out = random.sample(new_lines, random_sample)
            retained_lines = retained_lines + lines_out
        create_parent_directory(output_file)
        with open(output_file, "w") as f_out:
            f_out.write(f" {final_taxa_count} {str(get_site_count(input_file))}\n")
            for line in retained_lines:
                f_out.write(line)
