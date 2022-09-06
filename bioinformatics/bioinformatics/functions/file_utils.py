from logging import raiseExceptions
import os
import platform as pform
import shutil
import re
from typing import Optional
import subprocess


def list_files(input_folder: str) -> list[str]:
    file_list = []
    for root, dirs, files in os.walk(input_folder):
        for filename in files:
            # print(os.path.join(root, filename))
            file_list.append(os.path.join(root, filename))
    return file_list


def parse_output_file_name(in_file: str, output_folder: str = "output") -> str:
    suffix = suffix_parser(in_file)
    output_file = re.sub(rf"{suffix}$", f".out.{suffix}", in_file).replace(
        "..", "."
    )  # f".out{suffix}" > f".out.{suffix}"
    output_file = re.sub(r"/input/fasta/", f"/{output_folder}/", output_file)
    create_parent_directory(output_file)
    return output_file


# def execute_process(cmd: str) -> None:
#     subprocess.check_output(cmd).decode("utf-8")


def generate_process(
    cmd: str, stdout: Optional[bool] = False, output_result: Optional[bool] = False
) -> str:
    if stdout:
        subprocess.run(cmd)
        if output_result:
            return subprocess.check_output(cmd).decode("utf-8")
    else:
        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
        )
        if output_result:
            return subprocess.check_output(cmd).decode("utf-8")
    return None


def get_platform_info() -> dict[str, str]:
    bits, linkage = pform.architecture()
    info_dict = {
        "os": pform.system(),
        "platform": pform.platform(),
        "distribution": pform.freedesktop_os_release()["NAME"],
        "distribution_base": pform.freedesktop_os_release()["ID_LIKE"],
        "platform_release": pform.release(),
        "bits": bits,
        "linkage": linkage,
        "python_compiler": pform.python_compiler(),
        "python_implementation": pform.python_implementation(),
        "python_version": pform.python_version(),
    }
    if info_dict["os"] == "Linux":
        extra_info = generate_process(cmd="lscpu", output_result=True)
    elif info_dict["os"] == "Windows":
        extra_info = generate_process(cmd="systeminfo", output_result=True)
    else:
        extra_info = None
    for val in extra_info.split("\n"):
        parsed = val.split(":", 1)
        if len(parsed) == 2:
            info_dict[f"{parsed[0].strip()}"] = f"{parsed[1].strip()}"
    return info_dict


def create_parent_directory(directory: str) -> None:
    os.makedirs("/".join(directory.split("/")[0:-1]), exist_ok=True)


def resolve_suffix(suffix: str) -> str:
    if suffix in [".fa", ".fas", ".fasta", ".fna", ".fco", ".faa"]:
        suffix = "fasta"
    elif suffix in [".phy", ".phylip"]:
        suffix = "phylip"
    elif suffix in [".nex", ".nexus"]:
        suffix = "nexus"
    elif suffix in [".aln"]:
        suffix = "clustal"
    else:
        raise ValueError("Suffix Not Parsed")
    return suffix


def path_parser(in_file: str) -> str:
    split_path = os.path.splitext(in_file)
    current_path = split_path[0]
    current_suffix = split_path[1]
    return (current_path, current_suffix)


def suffix_parser(in_file: str, resolve: Optional[bool] = True) -> str:
    # suffix = re.search(r"[^.]+$", in_file).group()
    _, suffix = path_parser(in_file)
    if resolve:
        suffix = resolve_suffix(suffix)
    return suffix


def inject_parent_directory(in_file: str, injected_dir: str) -> str:
    injected_dir = re.sub(r"^\/|\/$", "", injected_dir)
    current_path, current_suffix = path_parser(in_file)
    current_parent = "/".join(current_path.split("/")[0:-1])
    current_file = current_path.split("/")[-1]
    new_parent = current_parent + f"/{injected_dir}/"
    new_file = f"{new_parent}{current_file}{current_suffix}"
    return new_file


def inject_grandparent_directory(in_file: str, injected_dir: str) -> str:
    injected_dir = re.sub(r"^\/|\/$", "", injected_dir)
    current_path, current_suffix = path_parser(in_file)
    current_parent = "/".join(current_path.split("/")[0:-2])
    current_file = current_path.split("/")[-1]
    new_parent = current_parent + f"/{injected_dir}/"
    new_file = f"{new_parent}{current_file}{current_suffix}"
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
    injected_prefix = re.sub(r"^\/|\/$", "", injected_prefix)
    injected_suffix = re.sub(r"^\/|\/$", "", injected_suffix)
    current_path, current_suffix = path_parser(in_file)
    current_parent = "/".join(current_path.split("/")[0:-1])
    current_file = current_path.split("/")[-1]
    new_file = injected_prefix + current_file + injected_suffix
    new_file = f"{current_parent}/{new_file}{current_suffix}"
    return new_file


def replace_suffix(in_file: str, new_suffix: str) -> str:
    suffix = suffix_parser(in_file, resolve=False)
    if new_suffix[0] != ".":
        new_suffix = "." + new_suffix
    new_file = in_file.replace(f"{suffix}", new_suffix)
    return new_file


def copy_directory(source: str, destination: str) -> None:
    shutil.copytree(source, destination)
