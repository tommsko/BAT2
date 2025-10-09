import logging
from enum import Enum

import MDAnalysis

from ...loader.simulation import Simulation, SegmentSupports
from .signature_generator import SignatureGenerator, SignatureGeneratorType, SignatureGenerationError
from ... import core_config

class NameGeneratorSource(Enum):
    CHAIN = 1
    MOLECULE = 2
    SEGMENT = 3

class NameGenerator(SignatureGenerator):
    def __init__(self, source: NameGeneratorSource) -> None:
        """
        Initializes NAME signature generator
        """
        if source == NameGeneratorSource.CHAIN:
            super().__init__("name_chain",
                             SignatureGeneratorType.DESCRIPTIVE,
                             "name")
        elif source == NameGeneratorSource.MOLECULE:
            super().__init__("name_molecule",
                             SignatureGeneratorType.DESCRIPTIVE,
                             "name")
        elif source == NameGeneratorSource.SEGMENT:
            super().__init__("name_segment",
                             SignatureGeneratorType.DESCRIPTIVE,
                             "name")
        else:
            raise RuntimeError("Unknown source for molecule names")

        self._source: NameGeneratorSource = source


    def _is_applicable(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> bool:
        """
        Checks whether the fragment can be described by its name
        :param simulation: simulation containing the fragment
        :param fragment: molecule to be described by its name
        :return: True if name signature is applicable, False otherwise
        """

        min_name_len: int = core_config.getint("signatures", "name_minimum_length")

        match self._source:
            case NameGeneratorSource.CHAIN:
                return (hasattr(fragment, "chainIDs")
                        and len(fragment.chainIDs) > 0
                        and len(fragment.chainIDs[0]) >= min_name_len)
            case NameGeneratorSource.MOLECULE:
                return (hasattr(fragment, "moltypes")
                        and len(fragment.moltypes) > 0
                        and len(fragment.moltypes[0]) >= min_name_len)
            case NameGeneratorSource.SEGMENT:
                return (hasattr(fragment, "segids")
                        and len(fragment.segids) > 0
                        and '_' in fragment.segids[0]
                        and len(fragment.segids[0].split("_", maxsplit=2)[-1]) >= min_name_len)
            case _:
                return False


    def _generate_signature(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> str:
        """
        Generates a name signature for the given fragment
        :param simulation: simulation containing the fragment
        :param fragment: molecule to be described by its name
        :raises SignatureGenerationError: if the generation fails
        :return: molecule name
        """

        match self._source:
            case NameGeneratorSource.CHAIN:
                return fragment.chainIDs[0]
            case NameGeneratorSource.MOLECULE:
                return fragment.moltypes[0]
            case NameGeneratorSource.SEGMENT:
                return fragment.segids[0].split("_", maxsplit=2)[-1]

        raise SignatureGenerationError("Unknown name source")


class ChainNameGenerator(NameGenerator):
    def __init__(self) -> None:
        """
        Initializes NAME signature generator w/ chain source (wrapper for NameGenerator)
        """
        super().__init__(source=NameGeneratorSource.CHAIN)

class MoleculeNameGenerator(NameGenerator):
    def __init__(self) -> None:
        """
        Initializes NAME signature generator w/ molecule source (wrapper for NameGenerator)
        """
        super().__init__(source=NameGeneratorSource.MOLECULE)

class SegmentNameGenerator(NameGenerator):
    def __init__(self) -> None:
        """
        Initializes NAME signature generator w/ segment source (wrapper for NameGenerator)
        """
        super().__init__(source=NameGeneratorSource.SEGMENT)