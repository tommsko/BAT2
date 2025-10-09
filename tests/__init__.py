import os.path
import shutil

TEST_DIR = "__test_delete_me"


def write_file(name: str, directory: str, content: str, append: bool = False) -> str:
    with open(os.path.join(directory, name), 'w' if not append else 'a') as fd:
        fd.write(content)
    return content


def check_file(name: str, directory: str, expected_content: str) -> bool:
    path: str = os.path.join(directory, name)
    assert os.path.exists(path), "target file does not exist"
    with open(path, 'r') as fd:
        content = fd.read()
        assert content == expected_content, f"target file content differs: got '{content}' expected '{expected_content}'"

    return True  # unnecessary, just for the consistency of tests


def delete_file(name: str, directory: str) -> None:
    path: str = os.path.join(directory, name)
    if os.path.exists(path):
        os.unlink(path)


def make_directory(name: str) -> None:
    if os.path.exists(name):
        return
    os.mkdir(name)


def delete_directory(name: str) -> None:
    shutil.rmtree(name)