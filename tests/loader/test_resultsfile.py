import json
import os.path

from .. import delete_file, check_file, write_file, delete_directory, make_directory, TEST_DIR
from ...src.loader.results_file import ResultsFile


# init empty
def test_resultsfile_init():
    make_directory(TEST_DIR)

    path = os.path.join(TEST_DIR, "test1.json")
    ResultsFile(path)

    assert os.path.exists(path)

    delete_directory(TEST_DIR)


# get simple
def test_resultsfile_get_simple():
    make_directory(TEST_DIR)

    path = os.path.join(TEST_DIR, "test2.json")
    write_file("test2.json", TEST_DIR, json.dumps({'hello': 1}))

    file = ResultsFile(path)
    assert file.get_item('hello') == 1

    delete_directory(TEST_DIR)


def test_load_empty():
    make_directory(TEST_DIR)

    path = os.path.join(TEST_DIR, "testX.json")
    write_file("testX.json", TEST_DIR, "")

    ResultsFile(path)  # this must not raise any exception

    delete_directory(TEST_DIR)


def test_resultsfile_exists():
    make_directory(TEST_DIR)

    path = os.path.join(TEST_DIR, "test2.json")
    write_file("test2.json", TEST_DIR, json.dumps({'hello': 1}))

    file = ResultsFile(path)
    assert file.item_exists("hello")
    assert not file.item_exists("world")

    delete_directory(TEST_DIR)


# get nested
def test_resultsfile_get_nested():
    make_directory(TEST_DIR)

    path = os.path.join(TEST_DIR, "test3.json")
    write_file("test3.json", TEST_DIR, json.dumps({'hello': {'world': 2}}))

    file = ResultsFile(path)
    assert file.get_item('hello/world') == 2

    delete_directory(TEST_DIR)


# set + get
def test_resultsfile_set_get():
    make_directory(TEST_DIR)

    path = os.path.join(TEST_DIR, "test4.json")
    file = ResultsFile(path)
    file.set_item('here/I/am', 42)

    file2 = ResultsFile(path)
    assert file2.get_item('here/I/am') == 42

    delete_directory(TEST_DIR)


# get default
def test_resultsfile_get_default():
    make_directory(TEST_DIR)

    path = os.path.join(TEST_DIR, "test5.json")
    file = ResultsFile(path)
    assert file.get_item('foo/bar/boo', default=68) == 68

    delete_directory(TEST_DIR)


# get type + incorrect
def test_resultsfile_get_type():
    make_directory(TEST_DIR)

    path = os.path.join(TEST_DIR, "test6.json")
    file = ResultsFile(path)
    file.set_item('foo/bar/boo', 68)

    assert file.get_item('foo/bar/boo') == 68
    assert file.get_item_safe('foo/bar/boo', int) == 68

    try:
        file.get_item_safe('foo/bar/boo', bool)
    except ValueError:
        ... # ok
    except Exception:
        assert False
    else:
        assert False

    delete_directory(TEST_DIR)