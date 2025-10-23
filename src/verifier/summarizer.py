import logging
import os.path
from collections import defaultdict
from enum import Enum

import MDAnalysis
from Bio.PDB import MMCIFIO, Structure, Model, Chain, Residue, Atom, PDBIO
from Bio.PDB.Polypeptide import is_aa

from ..loader.results_file import ResultsFile
from ..identifiers.identifier_manager import get_fragment_names
from ..loader.simulation import Simulation, SegmentSupports
from .verification_workers.verification_core import VerificationKind, Verifier
from .verification_workers.protein_nucleic_worker import ProteinNucleicVerificationWorker
from .verification_workers.compound_worker import CompoundVerificationWorker

logger: logging.Logger = logging.getLogger("verification")

class SegmentMoleculeType(Enum):
    PROTEIN = 1
    NUCLEIC = 2
    LIPID = 3
    CARBOHYDRATE = 4
    ATOM = 5
    WATER = 6
    UNKNOWN = 7

class ResultsSummarizer:
    # summarizes the results after the run of the pipeline

    @staticmethod
    def _classify_segment_type(sim: Simulation, results: ResultsFile, fragment_name: str) -> SegmentMoleculeType:
        """
        Classifies a specific segment (fragment) into one of the categories: protein, nucleic, lipid, carbohydrate, atom
        :param sim: simulation to export from
        :param results: results file associated with the simulation
        :param fragment_name: to be summarized
        :return: SegmentMoleculeType
        """
        is_protein = is_nucleic = is_lipid = is_carbohydrate = is_atom = is_water = is_unknown = False

        signatures: dict[str, str] = {}
        if results.item_exists(f"signatures/{fragment_name}/generated"):
            signatures = results.get_item(f"signatures/{fragment_name}/generated")

        if any('SMILES' in k.upper() and v.upper() == 'O' for k, v in signatures.items()):
            is_water = True
        elif sim.get_segment_flag(fragment_name, SegmentSupports.AMINO_ACID_RESIDUE_NAMES):
            is_protein = True
        elif sim.get_segment_flag(fragment_name, SegmentSupports.NUCLEIC_ACID_RESIDUE_NAMES):
            is_nucleic = True
        elif len(sim.get_fragment(fragment_name).atoms.ids) == 1:
            is_atom = True
        elif sim.get_segment_flag(fragment_name, SegmentSupports.ELEMENTS) or sim.fix_missing_elements(fragment_name):
            fragment: MDAnalysis.AtomGroup = sim.get_fragment(fragment_name)
            elements: list[str] = fragment.atoms.elements

            atom_c_cnt = sum(1 if e.upper() == "C" else 0 for e in elements)
            atom_o_cnt = sum(1 if e.upper() == "O" else 0 for e in elements)

            # TODO: this is very crude
            if atom_c_cnt > 10 and atom_o_cnt > 4 and 1 > (atom_c_cnt / atom_o_cnt) > 3:
                is_carbohydrate = True
            else:
                is_lipid = True
        else:
            is_unknown = True

        for field, value in [('protein', is_protein), ('nucleic', is_nucleic),
                             ('lipid', is_lipid), ('carbohydrate', is_carbohydrate),
                             ('atom', is_atom), ('water', is_water), ('unknown', is_unknown)]:
            results.set_item(f"summary/segments/{fragment_name}/macromolecule_type/{field}", value)

        if is_protein: return SegmentMoleculeType.PROTEIN
        if is_nucleic: return SegmentMoleculeType.NUCLEIC
        if is_lipid: return SegmentMoleculeType.LIPID
        if is_carbohydrate: return SegmentMoleculeType.CARBOHYDRATE
        if is_atom: return SegmentMoleculeType.ATOM
        if is_water: return SegmentMoleculeType.WATER
        return SegmentMoleculeType.UNKNOWN

    @staticmethod
    def _export_simulation_mmcif(results: ResultsFile, sim: Simulation, file_name: str = "simulation.mmcif",
                                 results_path: str = "summary/simulation_mmcif", filter_segment: str | None = None) -> None:
        """
        Exports the simulation into MMCIF file, inc. saving path
        :param results: results file associated with the simulation
        :param sim: to export
        :param file_name: to save mmcif as
        :param results_path: to update path of mmcif file in results
        :param filter_segment: to filter out segments (None for all segments)
        :return: None
        """

        def safe(val, default="."):
            if val is None or val.strip() == "":
                return default
            return val

        output_path: str = os.path.join(os.path.dirname(results.path), file_name)

        universe: MDAnalysis.Universe = sim.get_simulation()

        structure = Structure.Structure("0")
        model = Model.Model(0)
        structure.add(model)

        cached_residues: dict[int, dict[str, Residue.Residue]] = {}

        chains = {}
        for atom in universe.atoms:
            chain_id = atom.segid
            if filter_segment is not None and chain_id != filter_segment:
                continue

            resname = safe(atom.resname, default="UNK")
            atom_name = safe(atom.name, default="X")
            element = safe(atom.element, default="X")

            if chain_id not in chains:
                chain = Chain.Chain(chain_id)
                model.add(chain)
                chains[chain_id] = chain
            else:
                chain = chains[chain_id]

            if int(atom.resid) in cached_residues and resname in cached_residues[int(atom.resid)]:
                residue = cached_residues[int(atom.resid)][resname]
            else:
                residue = Residue.Residue((' ', int(atom.resid), ' '), resname, '')

                if str(atom.resid) not in cached_residues:
                    cached_residues[int(atom.resid)] = {}
                cached_residues[int(atom.resid)][resname] = residue
                chain.add(residue)

            bp_atom = Atom.Atom(atom_name, atom.position, 0.0, 0.0, '.', atom_name, atom.index, element=element.upper())
            residue.add(bp_atom)

        io = MMCIFIO()
        io.set_structure(structure)
        io.save(output_path)

        results.set_item(results_path, file_name)


    @staticmethod
    def summarize(sim: Simulation) -> Simulation:
        """
        Summarize the results of the verification.
        :param sim: simulation to be summarized
        :return: None
        """

        results: ResultsFile = sim.get_results_handler()

        # TODO: include this
        # if results.item_exists(f"summary"):
        #     return sim

        results.set_item(f"summary", {})

        try:
            ResultsSummarizer._export_simulation_mmcif(results, sim)
        except Exception as ex:
            logger.error(f"Unable to export simulation to mmcif file: {ex}")

        segments: list[str] = get_fragment_names(results)

        for segment in segments:

            try:
                ResultsSummarizer._export_simulation_mmcif(results, sim,
                                                           file_name=f"{segment}.mmcif",
                                                           results_path=f"summary/segments/{segment}/simulation_mmcif",
                                                           filter_segment=segment)
            except Exception as ex:
                logger.error(f"Unable to export segment '{segment}' to mmcif file: {ex}")

            segment_type: SegmentMoleculeType = SegmentMoleculeType.UNKNOWN
            try:
                segment_type = ResultsSummarizer._classify_segment_type(sim, results, segment)
            except Exception as ex:
                logger.error(f"Unable to classify type of segment '{segment}': {ex}")

            verifier: Verifier | None = None
            if segment_type in (SegmentMoleculeType.PROTEIN, SegmentMoleculeType.NUCLEIC):
                verifier = Verifier(sim, results, segment, VerificationKind.PROTEIN, ProteinNucleicVerificationWorker())
            elif segment_type in (SegmentMoleculeType.LIPID, SegmentMoleculeType.CARBOHYDRATE):
                verifier = Verifier(sim, results, segment, VerificationKind.PROTEIN, CompoundVerificationWorker())
            elif segment_type == SegmentMoleculeType.WATER:
                verifier = Verifier(sim, results, segment, VerificationKind.WATER, None)
            elif segment_type == SegmentMoleculeType.ATOM:
                verifier = Verifier(sim, results, segment, VerificationKind.ATOM, None)
            else:
                verifier = Verifier(sim, results, segment, VerificationKind.NONE, None)

            verified_identifiers: list[tuple[str, str, float, float | None]] = verifier.calculate_identifier_confidence()

            if verified_identifiers:
                results.set_item(f"summary/segments/{segment}/ident", verified_identifiers[0][0])
                results.set_item(f"summary/segments/{segment}/ident_type", verified_identifiers[0][1])
                results.set_item(f"summary/segments/{segment}/confidence", verified_identifiers[0][3])
                results.set_item(f"summary/segments/{segment}/verification", verified_identifiers)
            else:
                results.set_item(f"summary/segments/{segment}/ident", None)
                results.set_item(f"summary/segments/{segment}/ident_type", None)
                results.set_item(f"summary/segments/{segment}/confidence", None)
                results.set_item(f"summary/segments/{segment}/verification", None)

        return sim
