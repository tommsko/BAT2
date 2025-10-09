import hashlib
import logging
from enum import Enum

import MDAnalysis

from ...loader.simulation import Simulation, SegmentSupports
from .signature_generator import SignatureGenerator, SignatureGeneratorType, SignatureGenerationError
from ... import core_config


class FingerprintGenerator(SignatureGenerator):
    def __init__(self) -> None:
        """
        Initializes FINGERPRINT signature generator
        """
        super().__init__("fingerprint",
                         SignatureGeneratorType.FALLBACK,
                         "hash")


    def _is_applicable(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> bool:
        """
        Checks whether the fragment can be described by its fingerprint
        :param simulation: simulation containing the fragment
        :param fragment: molecule to be described by its fingerprint
        :return: True (always)
        """
        return True


    def _generate_signature(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> str:
        """
        Generates a fingerprint signature for the given fragment
        :param simulation: simulation containing the fragment
        :param fragment: molecule to be described by its fingerprint
        :raises SignatureGenerationError: if the generation fails
        :return: SHA512 hash
        """

        atom_count: int = fragment.n_atoms
        bonds: list[tuple[int, int]] = fragment.bonds.indices
        bonds.sort()
        return str(hashlib.sha512(f"{atom_count}_{bonds}".encode("utf-8")).hexdigest())

