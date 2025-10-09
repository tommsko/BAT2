import logging
import os.path

from .. import delete_file, check_file, write_file, delete_directory, make_directory, TEST_DIR
from ...src.loader.universe_loader import load, LoadStrategy, LOAD_STRATEGY_PATH
from ...src.loader.results_file import ResultsFile

def check_load(directory: str, results_file: ResultsFile):
    try:
        load(directory, results_file)
    except Exception as e:
        logging.error(f"load failed: {e}")
        return False
    else:
        return True

def test_itp_gro():
    make_directory(TEST_DIR)
    simulation_directory = "tests/data/16319"
    results_file = ResultsFile(os.path.join(TEST_DIR, "test1.json"))
    assert check_load(simulation_directory, results_file)
    assert results_file.get_item(LOAD_STRATEGY_PATH) == LoadStrategy.ITP_GRO.value
    delete_directory(TEST_DIR)

def test_top_gro():
    make_directory(TEST_DIR)
    simulation_directory = "tests/data/1010142"
    results_file = ResultsFile(os.path.join(TEST_DIR, "test2.json"))
    assert check_load(simulation_directory, results_file)
    assert results_file.get_item(LOAD_STRATEGY_PATH) in (LoadStrategy.TOP_GRO.value, LoadStrategy.ITP_GRO.value)
    delete_directory(TEST_DIR)

def test_tpr_gro():
    make_directory(TEST_DIR)
    simulation_directory = "tests/data/35193_a"
    results_file = ResultsFile(os.path.join(TEST_DIR, "test3.json"))
    assert check_load(simulation_directory, results_file)
    assert results_file.get_item(LOAD_STRATEGY_PATH) == LoadStrategy.TPR_GRO.value
    delete_directory(TEST_DIR)

def test_tpr():
    make_directory(TEST_DIR)
    simulation_directory = "tests/data/14976_a"
    results_file = ResultsFile(os.path.join(TEST_DIR, "test4.json"))
    assert check_load(simulation_directory, results_file)
    assert results_file.get_item(LOAD_STRATEGY_PATH) == LoadStrategy.TPR.value
    delete_directory(TEST_DIR)

def test_bad():
    make_directory(TEST_DIR)
    simulation_directory = "tests/data/14976"
    results_file = ResultsFile(os.path.join(TEST_DIR, "test5.json"))
    assert not check_load(simulation_directory, results_file)
    delete_directory(TEST_DIR)
