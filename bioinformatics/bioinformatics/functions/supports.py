import os
import re
from tqdm import tqdm

# input_folder = 'bioinformatics/input/test_files/supports'


def rearrange_supports(input_folder: str, **kwargs) -> None:
    for root, dirs, files in os.walk(input_folder):
        for filename in tqdm(files):
            rearrange_support(os.path.join(root, filename), **kwargs)


def rearrange_support(
    in_file: str, output_suffix: str, bs_order: list[int], boot_tree_sep: str = ";"
) -> None:
    # in_file = 'bioinformatics/input/test_files/supports/L1.treefile'
    # in_file = 'bioinformatics/input/test_files/supports/L1.ufboot'
    if len(bs_order) not in [0, 1]:
        search_string = "\)" + "([0-9|\.]+)\/" * (len(bs_order) - 1) + "([0-9|\.]+)\:"
        replace_string = ")\\" + "/\\".join([str(x) for x in bs_order]) + ":"
        filename = in_file.split(".")[0]
        suffix = in_file.split(".")[-1]
        with open(in_file, "r") as f:  # read tree
            contents = f.read()
        output = []
        if suffix in ["ufboot"]:
            trees = contents.split(boot_tree_sep)
        elif suffix in ["tre", "treefile"]:
            trees = [contents]
        for tree in trees:
            if tree not in ["", "\n"]:
                new_tree = re.sub(
                    rf"{search_string}",
                    rf"{replace_string}",
                    tree,
                )
                output.append(new_tree)
        output = f"{boot_tree_sep}".join(output)
        with open(filename + output_suffix + "." + suffix, "w") as f:  # write results
            f.write(output)
    else:
        raise ValueError("2 or more branch supports per branch required.")
