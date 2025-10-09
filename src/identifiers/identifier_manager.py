import logging
from enum import Enum
from logging import Logger

from .. import core_config
from ..loader.results_file import ResultsFile
from ..loader.simulation import Simulation
from ..signatures.signature_manager import SignatureType

from .identifier_resolvers.identifier_resolver import IdentifierResolver, IdentifierResolutionError
from .identifier_resolvers.blast import BLASTResolverPDBAA, BLASTResolverPDBNT, BLASTResolverUNIPROT
from .identifier_resolvers.pdb_struct import PDBStructResolverProtein, PDBStructResolverNucleotide
from .identifier_resolvers.alphafind_struct import AlphaFindStructResolver
from .identifier_resolvers.chembl import CHEMBLResolver
from .identifier_resolvers.pubchem import PUBCHEMResolver
from .identifier_resolvers.atom import AtomResolver
from .identifier_resolvers.name import NameResolverCCD, NameResolverMAPS, NameResolverMANUAL
from .identifier_resolvers.idsm import IDSMResolverPubchem, IDSMResolverChembl, IDSMResolverChebi

logger: Logger = logging.getLogger("identifier")

class IdentifierType(Enum):  # order here represents order in which identifiers are attempted
    BLAST_PDB_PROT = BLASTResolverPDBAA
    BLAST_PDB_NT = BLASTResolverPDBNT
    BLAST_UNIPROT = BLASTResolverUNIPROT
    PDB_STRUCT_PROTEIN = PDBStructResolverProtein
    PDB_STRUCT_NUCLEOTIDE = PDBStructResolverNucleotide
    ALPHAFIND_STRUCT = AlphaFindStructResolver
    IDSM_PUBCHEM = IDSMResolverPubchem
    IDSM_CHEMBL = IDSMResolverChembl
    IDSM_CHEBI = IDSMResolverChebi
    CHEMBL_STRUCT = CHEMBLResolver
    PUBCHEM_STRUCT = PUBCHEMResolver
    ATOM = AtomResolver
    NAME_CCD = NameResolverCCD
    NAME_MAPS = NameResolverMAPS
    NAME_MANUAL = NameResolverMANUAL


def map_signature_generator_to_signature_type() -> dict[str, str]:
    """
    :return: associated name of signature generator (found in results file) to its signature type
    """
    return {
        generator.name: (generator.value)().signature_type for generator in SignatureType
    }


def get_fragment_names(results: ResultsFile) -> list[str]:
    """
    Extracts fragments stored inside the results file (from signature generation step)
    :param results: containing fragments
    :return: list of fragment names
    """
    return list(results.get_item("signatures").keys())


def get_generated_signatures(results: ResultsFile, fragment_name: str) -> dict[str, str]:
    """
    Extracts signatures stored in results file
    :param results: containing signatures
    :param fragment_name: to extract signatures of
    :return: mapping of signature generator to signature string
    """
    return results.get_item(f"signatures/{fragment_name}/generated").items()


def check_init_identifier_data(results: ResultsFile, fragment_name: str, signature_generator_name: str) -> None:
    """
    Checks whether identifier data segment in results file is ready, if not creates empty one
    :param results: containing signatures/identifier
    :param fragment_name: for which the data field is concerned
    :param signature_generator_name: for which data field is concerned
    :return: None (updates results file)
    """
    if not results.item_exists(f"identifiers/{fragment_name}/{signature_generator_name}/attempted"):
        results.set_item(f"identifiers/{fragment_name}/{signature_generator_name}/attempted", [])


def is_signature_identified(results: ResultsFile, fragment_name: str, signature_generator_name: str) -> bool:
    """
    Checks whether a given signature is identified on either similarity or identity level
    :param results: containing signatures/identifier
    :param fragment_name: for which the data field is concerned
    :param signature_generator_name: for which data field is concerned
    :return: True if signature is identified on either identity or similarity level
    """

    if not results.item_exists(f"identifiers/{fragment_name}/{signature_generator_name}/"):
        return False

    possible_identifiers = list(results.get_item(f"identifiers/{fragment_name}/{signature_generator_name}/").keys())
    possible_identifiers = [k for k in possible_identifiers if k != "attempted"]

    for identifier in possible_identifiers:
        if not results.item_exists(f"identifiers/{fragment_name}/{signature_generator_name}/{identifier}/hits"):
            return False
        if len(results.get_item(f"identifiers/{fragment_name}/{signature_generator_name}/{identifier}/hits")) > 0:
            return True

    return False


def add_attempted_identifier(results: ResultsFile, fragment_name: str, signature_generator_name: str,
                             resolver_name: str) -> None:
    """
    Saves attempted resolver into results file
    :param results: to store attempted identifier resolver
    :param fragment_name: for which resolution was attempted
    :param signature_generator_name: for which resolution was attempted
    :param resolver_name: name of identifier resolver that was attempted
    :return: None (updates results file)
    """
    curr: list[str] = results.get_item(f"identifiers/{fragment_name}/{signature_generator_name}/attempted")
    if resolver_name in curr:
        return

    curr.append(resolver_name)
    results.set_item(f"identifiers/{fragment_name}/{signature_generator_name}/attempted", curr)

def add_identifiers(results: ResultsFile, fragment_name: str, signature_generator_name: str, signature: str,
                    resolver_name: str, identifier_result_type: str, identifiers: list[str] | list[tuple[str, float]], is_similarity: bool) -> None:
    """
    Saves successful identifier resolution into results file
    :param results: to store resolved identifier
    :param fragment_name: for which resolution was done
    :param signature_generator_name: for which resolution was done
    :param resolver_name: of which succeeded in identification
    :param identifier_result_type: type of identifier hits
    :param identifiers: associated identifiers to fragment/signature
    :param is_similarity: if True, saves as similarity identification opposed to exact-match
    :return: None (updates results file)
    """
    results.set_item(f"identifiers/{fragment_name}/{signature_generator_name}/{resolver_name}", {
        "signature": signature,
        "hit_type": identifier_result_type,
        "hits": []
    })

    if is_similarity:
        identifiers = sorted(identifiers, key=lambda tup: tup[1], reverse=True)
        identifiers = [(k, round(float(v) * 100, 2) if float(v) <= 1.0 else v) for k, v in identifiers]
        results.set_item(f"identifiers/{fragment_name}/{signature_generator_name}/{resolver_name}/hits", identifiers)
    else:
        identifiers = [(k, 100.0) for k in identifiers]
        results.set_item(f"identifiers/{fragment_name}/{signature_generator_name}/{resolver_name}/hits", identifiers)


class IdentifierResolutionManager:

    generator_name_to_signature_type: dict[str, str] = map_signature_generator_to_signature_type()

    def __init__(self) -> None:
        """
        Initializes the IdentifierResolutionManager
        """

        self._overwrite_existing: bool = core_config.getboolean('identifiers', 'overwrite_existing')

    @staticmethod
    def resolve_signature(identifier: IdentifierResolver, identifier_name: str, identifier_result_type: str,
                          signature: str, signature_generator_name: str,
                          results: ResultsFile, segment: str,
                          is_similarity: bool) -> bool:
        """
        Attempts to resolve signature given
        :param identifier: to use for resolution
        :param identifier_name: name of the above
        :param identifier_result_type: type of identifier generated
        :param signature: to be resolved
        :param signature_generator_name: name of the above
        :param results: to store identifiers and metadata
        :param segment: of which identification is being done
        :param is_similarity: if True, similarity identification is attempted instead of exact-match
        :return: True if identification was successful
        """
        if identifier.is_applicable(signature,
                                    IdentifierResolutionManager.
                                            generator_name_to_signature_type[signature_generator_name],
                                    as_similarity=is_similarity):
            add_attempted_identifier(results, segment, signature_generator_name, identifier_name)
        else:
            return False

        if (idents := identifier.resolve_signature(signature,
                                                   as_similarity=is_similarity,
                                                   working_directory=results.working_directory)) is not None:
            add_identifiers(results, segment, signature_generator_name, signature,
                            identifier_name, identifier_result_type, idents, is_similarity=is_similarity)
            return True

        return False

    def resolve_identifiers(self, simulation: Simulation) -> Simulation:
        """
        Attempts to resolve identifiers for all signatures in the simulation
        Beware, fragments are extracted from results file, not the simulation
        :param results_file_path: path to the results file of the simulation
        :param simulation: (unused, only for task chaining)
        :return: the same simulation (for task chaining)
        """

        logger.info("Begin identification of signatures for simulation '%s'" % simulation.simulation_directory)

        results: ResultsFile = simulation.get_results_handler()
        segments: list[str] = get_fragment_names(results)

        for segment in segments:
            logger.info(f"Begin identification of signatures for simulation '{simulation.simulation_directory}', "
                        f"segment '{segment}'")

            for signature_generator_name, signature in get_generated_signatures(results, segment):
                check_init_identifier_data(results, segment, signature_generator_name)

                if not self._overwrite_existing and is_signature_identified(results, segment, signature_generator_name):
                    continue

                for identifier_type in IdentifierType:
                    identifier_name: str = identifier_type.name
                    identifier: IdentifierResolver = (identifier_type.value)()
                    identifier_result_type: str = identifier.identifier_type

                    logger.debug(f"Attempting to identify signatures from '{signature_generator_name}' "
                                 f"using '{identifier_name}'")

                    # identity identification
                    if IdentifierResolutionManager.resolve_signature(identifier, identifier_name, identifier_result_type,
                                                                     signature, signature_generator_name,
                                                                     results, segment, is_similarity=False):
                        logger.debug("Successful identity identification!")
                        #break # we want everything

                else:  # no exact-match identifier found, continue with similarity
                    for identifier_type in IdentifierType:
                        identifier_name: str = identifier_type.name
                        identifier: IdentifierResolver = (identifier_type.value)()
                        identifier_result_type: str = identifier.identifier_type

                        # similarity identification
                        if IdentifierResolutionManager.resolve_signature(identifier, identifier_name, identifier_result_type,
                                                                         signature, signature_generator_name,
                                                                         results, segment, is_similarity=True):
                            logger.debug("Successful similarity identification!")
                            #break we want everything
                    #logger.debug("Failed identification!")

        return simulation

