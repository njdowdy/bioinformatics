import os
from typing import Optional


def create_parent_directory(directory: str) -> None:
    os.makedirs("/".join(directory.split("/")[0:-1]), exist_ok=True)


def suffix_parser(in_file: str) -> str:
    # suffix = re.search(r"[^.]+$", in_file).group()
    suffix = os.path.splitext(in_file)[1]
    return suffix


def inject_parent_directory(in_file: str, injected_dir: str) -> str:
    injected_dir = injected_dir.replace("/", "")
    current_path = os.path.splitext(in_file)[0]
    current_parent = "/".join(current_path.split("/")[0:-1])
    current_file = current_path.split("/")[-1]
    new_parent = current_parent + f"/{injected_dir}/"
    new_file = f"{new_parent}{current_file}"
    return new_file


def replace_parent_directory(in_file: str, search_dir: str, replace_dir: str) -> str:
    new_file = in_file.replace(f"{search_dir}", f"{replace_dir}")
    return new_file


def combine_files(directory: str, out_file: str) -> None:
    create_parent_directory(out_file)
    for root, dirs, files in os.walk(directory):
        with open(out_file, "w") as outfile:
            for filename in files:
                with open(os.path.join(root, filename)) as infile:
                    for line in infile:
                        outfile.write(line)
    return None


def inject_prefix_suffix(
    in_file: str,
    injected_prefix: Optional[str] = "",
    injected_suffix: Optional[str] = "",
) -> str:
    injected_prefix = injected_prefix.replace("/", "")
    injected_suffix = injected_suffix.replace("/", "")
    current_path = os.path.splitext(in_file)[0]
    current_parent = "/".join(current_path.split("/")[0:-1])
    current_file = current_path.split("/")[-1]
    new_file = injected_prefix + current_file + injected_suffix
    new_file = f"{current_parent}{new_file}"
    return new_file
