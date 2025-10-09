import os.path
import string
from os import listdir
from os.path import isfile, join, isdir
import random

from rich import print


def list_files_in_directory(path: str, only_basename: bool = False) -> list[str]:
    """
    Lists files in the specified directory
    :param path: to the directory
    :param only_basename: if true, only files' basename are returned, otherwise full absolute paths
    :return: list of files
    """

    return [file if only_basename else
            os.path.abspath(os.path.join(path, file))
            for file in listdir(path) if isfile(join(path, file))]


def filter_files_by_extension(files: list[str], allowed_extensions: list[str]) -> list[str]:
    """
    Filters files ending with allowed extensions from all files given
    :param files: list of paths
    :param allowed_extensions: list of allowed extensions, case-insensitive
    :return: list of files with allowed extensions
    """
    allowed_extensions = [x.lower() for x in allowed_extensions]

    matching_files: list[str] = []
    for file in files:
        extension = file.lower().split('.')[-1]
        if extension in allowed_extensions:
            matching_files.append(file)
    return matching_files


def split_files_by_extension(files: list[str], extension: str) -> tuple[list[str], list[str]]:
    """
    Splits files into two groups, those ending with specific extension and those which do not
    :param files: list of paths
    :param extension: to split by
    :return: list of files with matching extension and those without
    """

    files_matching_extension: list[str] = filter_files_by_extension(files, [extension])
    return files_matching_extension, list(set(files) - set(files_matching_extension))

def list_folders_in_directory(path: str, stdout_progress: bool = True) -> list[str]:
    """
    Traverses directory and returns all folders found
    :param path: of the directory
    :param stdout_progress: if true, print progress to stdout
    :return: list of absolute paths of folders
    """

    if stdout_progress:
        print("[cyan][b]-> Listing simulations to process...[/b][/cyan]", end="\r")

    ret: list[str] = [
        os.path.abspath(join(path, subpath))
        for subpath in listdir(path)
        if isdir(join(path, subpath))
    ]
    if stdout_progress:
        print(f"[green][b]✔  Found {len(ret)} simulations              [/b][/green]")

    return ret

def mk_random_filename(extension: str, dir_path: str, name_len: int = 30) -> str:
    """
    Generates a random filename with given extension.
    :param extension: of the file, without the dot
    :param dir_path: directory where the file should be saved
    :param name_len: length of the name (excluding extension)
    :return: path to such file
    """
    extension = extension.replace(".", "")  # because of course we forget to

    # https://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits
    filename: str = "".join(
        random.choice(string.ascii_letters + string.digits) for _ in range(name_len)
    )

    return os.path.abspath(os.path.join(dir_path, f"{filename}.{extension}"))