import logging
from enum import Enum

from .. import core_config
from ..loader.simulation import Simulation
from ..loader.results_file import ResultsFile
from .signature_generators.signature_generator import SignatureGenerator
from .signature_generators.smiles import SimpleSmilesGenerator, SimpleSmilesGeneratorNoSanitize
from .signature_generators.atom import AtomGenerator
from .signature_generators.fasta import MapFastaGeneratorProtein, MapFastaGeneratorNucleotide, ResnameFastaGeneratorProtein, ResnameFastaGeneratorNucleotide
from .signature_generators.pdb import (PDBGeneratorSanitizeProtein, PDBGeneratorNoSanitizeProtein, PDBBackboneGeneratorSanitizeProtein, PDBBackboneGeneratorNoSanitizeProtein,
                                       PDBGeneratorSanitizeNucleic, PDBGeneratorNoSanitizeNucleic, PDBBackboneGeneratorSanitizeNucleic, PDBBackboneGeneratorNoSanitizeNucleic)
from .signature_generators.name import ChainNameGenerator, MoleculeNameGenerator, SegmentNameGenerator
from .signature_generators.fingerprint import FingerprintGenerator

logger: logging.Logger = logging.getLogger("signature")

class SignatureType(Enum):  # order here represents order in which signatures are generated
    FASTA_MAP_NUCLEIC = MapFastaGeneratorNucleotide
    FASTA_RESNAME_NUCLEIC = ResnameFastaGeneratorNucleotide
    FASTA_MAP_PROTEIN = MapFastaGeneratorProtein
    FASTA_RESNAME_PROTEIN = ResnameFastaGeneratorProtein
    ATOM = AtomGenerator
    SMILES = SimpleSmilesGenerator
    SMILES_NOT_SANITIZED = SimpleSmilesGeneratorNoSanitize
    PDB_SANITIZED_PROTEIN = PDBGeneratorSanitizeProtein
    PDB_BACKBONE_SANITIZED_PROTEIN = PDBBackboneGeneratorSanitizeProtein
    PDB_NOT_SANITIZED_PROTEIN = PDBGeneratorNoSanitizeProtein
    PDB_BACKBONE_NOT_SANITIZED_PROTEIN = PDBBackboneGeneratorNoSanitizeProtein
    PDB_SANITIZED_NUCLEIC = PDBGeneratorSanitizeNucleic
    PDB_BACKBONE_SANITIZED_NUCLEIC = PDBBackboneGeneratorSanitizeNucleic
    PDB_NOT_SANITIZED_NUCLEIC = PDBGeneratorNoSanitizeNucleic
    PDB_BACKBONE_NOT_SANITIZED_NUCLEIC = PDBBackboneGeneratorNoSanitizeNucleic
    NAME_CHAIN = ChainNameGenerator
    NAME_MOLECULE = MoleculeNameGenerator
    NAME_SEGMENT = SegmentNameGenerator
    FINGERPRINT = FingerprintGenerator

class SignatureGeneratorManager:
    def __init__(self) -> None:
        """
        Initializes the SignatureGeneratorManager
        """

        self._overwrite_existing: bool = core_config.getboolean('signatures', 'overwrite_existing')
        self._signature_gen_attempts: int = core_config.getint('signatures', 'generation_attempts')

    def generate_signatures(self, simulation: Simulation) -> Simulation:
        """
        Attempts to generate signatures for all fragments in the simulation
        :param simulation: to describe fragments of
        :return: The same simulation (for task chaining)
        """

        logger.info("Begin generating signatures for simulation '%s'" % simulation.simulation_directory)

        results_file: ResultsFile = simulation.get_results_handler()

        segments: set[str] = simulation.get_segments()
        for segment in segments:
            if self._overwrite_existing or not results_file.item_exists(f"signatures/{segment}"):
                results_file.set_item(f"signatures/{segment}/attempted", [])
                results_file.set_item(f"signatures/{segment}/generated", {})

        for segment_name in segments:
            logger.info(f"Begin generating signatures for simulation '{simulation.simulation_directory}', "
                        f"segment '{segment_name}'")

            for signature_type in SignatureType:

                if (not self._overwrite_existing
                        and signature_type.name in results_file.get_item(f"signatures/{segment_name}/generated")):
                    continue

                logger.debug(f"Attempting to generate signatures for segment '{segment_name}' "
                             f"using '{signature_type.name}'")

                generator: SignatureGenerator = (signature_type.value)()
                if generator.is_applicable(simulation, simulation.get_fragment(segment_name)):
                    current_attempted: list[str] = results_file.get_item(f"signatures/{segment_name}/attempted")
                    results_file.set_item(f"signatures/{segment_name}/attempted",
                                          current_attempted + [signature_type.name])
                else:
                    logger.debug(f"... '{signature_type.name}' is not applicable")
                    continue

                for _ in range(self._signature_gen_attempts):
                    if ((signature := generator.generate_signature(simulation, simulation.get_fragment(segment_name)))
                            is not None):
                        current_signs: dict[str, str] = results_file.get_item(f"signatures/{segment_name}/generated")
                        current_signs[signature_type.name] = signature
                        results_file.set_item(f"signatures/{segment_name}/generated", current_signs)
                        logger.debug(f"... '{signature_type.name}' generated signature {signature[:10]}...")
                        break
                else:
                    logger.debug(f"... '{signature_type.name}' failed to generate signature")

        return simulation
