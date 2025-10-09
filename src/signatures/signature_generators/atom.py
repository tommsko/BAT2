import logging

import MDAnalysis

from ...loader.simulation import Simulation, SegmentSupports
from .signature_generator import SignatureGenerator, SignatureGeneratorType, SignatureGenerationError
from ...constants import ATOMTABLE_ELEMENTS_MASSES


class AtomGenerator(SignatureGenerator):
    def __init__(self) -> None:
        """
        Initializes ATOM signature generator
        """
        super().__init__("atom",
                         SignatureGeneratorType.ANALYTICAL,
                         "element")



    def _is_applicable(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> bool:
        """
        Checks whether the fragment can be described by a single atomic element
        :param simulation: simulation containing the fragment
        :param fragment: ~~molecule~~ atom to be described by a single atomic element
        :return: True if atom signature is applicable, False otherwise
        """
        fragment_name: str = fragment.segids[0]

        # we need atomic elements
        if not simulation.get_segment_flag(fragment_name, SegmentSupports.ELEMENTS):
            if not simulation.fix_missing_elements(fragment_name):
                return False

        # we need exactly one atom
        return len(fragment.atoms.ids) == 1


    def _generate_signature(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> str:
        """
        Generates an atomic element signature for the given fragment
        :param simulation: simulation containing the fragment
        :param fragment: ~~molecule~~ atom to be described as a single atomic element
        :raises SignatureGenerationError: if the generation fails
        :return: atomic element (Title-case)
        """

        external_elements: list[str] | None = None
        if any(x == '' for x in fragment.atoms.elements):
            # non-deterministic bug in MDAnalysis when newly set elements 'disappear'
            external_elements = []
            simulation.fix_missing_elements(fragment.segids[0], external_elements)

        element: str = fragment.atoms.elements[0] if external_elements is None else external_elements[0]
        element = element.title()

        # sanitize
        if element not in ATOMTABLE_ELEMENTS_MASSES.keys():
            raise SignatureGenerationError(f"Invalid atomic element: {element}")

        return element

