import os
from .. import delete_file, check_file, write_file, delete_directory, make_directory, TEST_DIR
from ...src.loader.simulation import Simulation, SegmentSupports, SimulationParameter


def test_sim_14976_a():
    simulation_directory = "tests/data/14976_a"

    make_directory(TEST_DIR)
    simulation = Simulation(simulation_directory, os.path.join(TEST_DIR, "test1.json"))
    simulation.load_simulation()

    assert simulation.get_simulation_parameter(SimulationParameter.ATOM_COUNT) == 82308
    assert simulation.get_simulation_parameter(SimulationParameter.RESIDUE_COUNT) == 14896
    assert simulation.get_simulation_parameter(SimulationParameter.SEGMENT_COUNT) == 4
    assert simulation.get_simulation_parameter(SimulationParameter.FRAGMENT_COUNT) == 14896

    assert simulation.get_segments() == {'seg_0_POPC', 'seg_2_NA', 'seg_1_SOL', 'seg_3_CL'}

    assert simulation.get_segment_flag("seg_0_POPC", SegmentSupports.ELEMENTS) is False
    assert simulation.get_segment_flag("seg_0_POPC", SegmentSupports.MASSES) is True
    assert simulation.get_segment_flag("seg_0_POPC", SegmentSupports.POSITIONS) is False
    assert simulation.get_segment_flag("seg_0_POPC", SegmentSupports.NAMES) is True

    assert simulation.get_segment_flag("seg_2_NA", SegmentSupports.ELEMENTS) is True
    assert simulation.get_segment_flag("seg_2_NA", SegmentSupports.MASSES) is True
    assert simulation.get_segment_flag("seg_2_NA", SegmentSupports.POSITIONS) is False
    assert simulation.get_segment_flag("seg_2_NA", SegmentSupports.NAMES) is True

    assert simulation.get_segment_flag("seg_1_SOL", SegmentSupports.ELEMENTS) is True
    assert simulation.get_segment_flag("seg_1_SOL", SegmentSupports.MASSES) is True
    assert simulation.get_segment_flag("seg_1_SOL", SegmentSupports.POSITIONS) is False
    assert simulation.get_segment_flag("seg_1_SOL", SegmentSupports.NAMES) is True

    assert simulation.get_segment_flag("seg_3_CL", SegmentSupports.ELEMENTS) is True
    assert simulation.get_segment_flag("seg_3_CL", SegmentSupports.MASSES) is True
    assert simulation.get_segment_flag("seg_3_CL", SegmentSupports.POSITIONS) is False
    assert simulation.get_segment_flag("seg_3_CL", SegmentSupports.NAMES) is True

    delete_directory(TEST_DIR)

def test_sim_16319():
    simulation_directory = "tests/data/16319"

    make_directory(TEST_DIR)
    simulation = Simulation(simulation_directory, os.path.join(TEST_DIR, "test1.json"))
    simulation.load_simulation()

    assert simulation.get_simulation_parameter(SimulationParameter.ATOM_COUNT) == 12208
    assert simulation.get_simulation_parameter(SimulationParameter.RESIDUE_COUNT) == 2952
    assert simulation.get_simulation_parameter(SimulationParameter.SEGMENT_COUNT) == 2952
    assert simulation.get_simulation_parameter(SimulationParameter.FRAGMENT_COUNT) == 2952

    assert simulation.get_segments() == {'DPPC', 'SOL', 'Na', 'Cl'}

    assert simulation.get_segment_flag("DPPC", SegmentSupports.ELEMENTS) is True
    assert simulation.get_segment_flag("DPPC", SegmentSupports.MASSES) is True
    assert simulation.get_segment_flag("DPPC", SegmentSupports.POSITIONS) is True
    assert simulation.get_segment_flag("DPPC", SegmentSupports.NAMES) is True

    assert simulation.get_segment_flag("SOL", SegmentSupports.ELEMENTS) is True
    assert simulation.get_segment_flag("SOL", SegmentSupports.MASSES) is True
    assert simulation.get_segment_flag("SOL", SegmentSupports.POSITIONS) is True
    assert simulation.get_segment_flag("SOL", SegmentSupports.NAMES) is True

    assert simulation.get_segment_flag("Na", SegmentSupports.ELEMENTS) is True
    assert simulation.get_segment_flag("Na", SegmentSupports.MASSES) is True
    assert simulation.get_segment_flag("Na", SegmentSupports.POSITIONS) is True
    assert simulation.get_segment_flag("Na", SegmentSupports.NAMES) is True

    assert simulation.get_segment_flag("Cl", SegmentSupports.ELEMENTS) is True
    assert simulation.get_segment_flag("Cl", SegmentSupports.MASSES) is True
    assert simulation.get_segment_flag("Cl", SegmentSupports.POSITIONS) is True
    assert simulation.get_segment_flag("Cl", SegmentSupports.NAMES) is True

    delete_directory(TEST_DIR)

def test_sim_35193_a():
    simulation_directory = "tests/data/35193_a"

    make_directory(TEST_DIR)
    simulation = Simulation(simulation_directory, os.path.join(TEST_DIR, "test1.json"))
    simulation.load_simulation()

    assert simulation.get_simulation_parameter(SimulationParameter.ATOM_COUNT) == 53962
    assert simulation.get_simulation_parameter(SimulationParameter.RESIDUE_COUNT) == 9362
    assert simulation.get_simulation_parameter(SimulationParameter.SEGMENT_COUNT) == 4
    assert simulation.get_simulation_parameter(SimulationParameter.FRAGMENT_COUNT) == 9362

    assert simulation.get_segments() == {'seg_2_Na_s', 'seg_1_TIP3p', 'seg_0_POPC', 'seg_3_Cl_s'}

    assert simulation.get_segment_flag("seg_2_Na_s", SegmentSupports.ELEMENTS) is True
    assert simulation.get_segment_flag("seg_2_Na_s", SegmentSupports.MASSES) is True
    assert simulation.get_segment_flag("seg_2_Na_s", SegmentSupports.POSITIONS) is True
    assert simulation.get_segment_flag("seg_2_Na_s", SegmentSupports.NAMES) is True

    assert simulation.get_segment_flag("seg_1_TIP3p", SegmentSupports.ELEMENTS) is True
    assert simulation.get_segment_flag("seg_1_TIP3p", SegmentSupports.MASSES) is True
    assert simulation.get_segment_flag("seg_1_TIP3p", SegmentSupports.POSITIONS) is True
    assert simulation.get_segment_flag("seg_1_TIP3p", SegmentSupports.NAMES) is True

    assert simulation.get_segment_flag("seg_0_POPC", SegmentSupports.ELEMENTS) is True
    assert simulation.get_segment_flag("seg_0_POPC", SegmentSupports.MASSES) is True
    assert simulation.get_segment_flag("seg_0_POPC", SegmentSupports.POSITIONS) is True
    assert simulation.get_segment_flag("seg_0_POPC", SegmentSupports.NAMES) is True

    assert simulation.get_segment_flag("seg_3_Cl_s", SegmentSupports.ELEMENTS) is True
    assert simulation.get_segment_flag("seg_3_Cl_s", SegmentSupports.MASSES) is True
    assert simulation.get_segment_flag("seg_3_Cl_s", SegmentSupports.POSITIONS) is True
    assert simulation.get_segment_flag("seg_3_Cl_s", SegmentSupports.NAMES) is True

    delete_directory(TEST_DIR)

def test_sim_1010142():
    simulation_directory = "tests/data/1010142"

    make_directory(TEST_DIR)
    simulation = Simulation(simulation_directory, os.path.join(TEST_DIR, "test1.json"))
    simulation.load_simulation()

    assert simulation.get_simulation_parameter(SimulationParameter.ATOM_COUNT) == 27755
    assert simulation.get_simulation_parameter(SimulationParameter.RESIDUE_COUNT) == 8860
    assert simulation.get_simulation_parameter(SimulationParameter.SEGMENT_COUNT) == 8769
    assert simulation.get_simulation_parameter(SimulationParameter.FRAGMENT_COUNT) == 8769

    assert simulation.get_segments() == {'SOL', 'Protein'}

    assert simulation.get_segment_flag("SOL", SegmentSupports.ELEMENTS) is True
    assert simulation.get_segment_flag("SOL", SegmentSupports.MASSES) is True
    assert simulation.get_segment_flag("SOL", SegmentSupports.POSITIONS) is True
    assert simulation.get_segment_flag("SOL", SegmentSupports.NAMES) is True

    assert simulation.get_segment_flag("Protein", SegmentSupports.ELEMENTS) is True
    assert simulation.get_segment_flag("Protein", SegmentSupports.MASSES) is True
    assert simulation.get_segment_flag("Protein", SegmentSupports.POSITIONS) is True
    assert simulation.get_segment_flag("Protein", SegmentSupports.NAMES) is True

    delete_directory(TEST_DIR)