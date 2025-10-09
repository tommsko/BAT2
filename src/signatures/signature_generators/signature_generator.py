import logging
from enum import Enum
from logging import Logger

from ... import core_config
from ...loader.simulation import Simulation, SegmentSupports
import MDAnalysis

logger: Logger = logging.getLogger("signature")

class SignatureGenerationError(Exception):
    def __init__(self, message):
        super().__init__(message)


class SignatureGeneratorType(Enum):
    ANALYTICAL = 0
    DESCRIPTIVE = 1
    FALLBACK = 2


class SignatureGenerator:
    """
    Abstract class for signature generators which input MDAnalysis fragment and output a signature (string)
    """

    def __init__(self, generator_name: str, generator_type: SignatureGeneratorType, signature_type: str) -> None:
        """
        Initializes the SignatureGenerator with the given signature type
        :param generator_name: string name of the signature generator
        :param generator_type: type of the signature generator
        :param signature_type: string type of the signature
        """
        self.generator_name: str = generator_name
        self.generator_type: SignatureGeneratorType = generator_type
        self.signature_type: str = signature_type

        self._is_allowed: bool = core_config.getboolean('signatures', generator_name)

    def is_applicable(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> bool:
        """
        Checks whether the signature generator is applicable to the given fragment
        :param simulation: simulation containing the fragment
        :param fragment: molecule to be described by the signature generator
        :return: True if signature generator is applicable, False otherwise
        """
        if not self._is_allowed:
            return False
        return self._is_applicable(simulation, fragment)

    def _is_applicable(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> bool:
        """
        Implementation of check whether the signature generator is applicable to the given fragment
        :param simulation: simulation containing the fragment
        :param fragment: molecule to be described by the signature generator
        :return: True if signature generator is applicable, False otherwise
        """
        raise NotImplementedError

    def generate_signature(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> str | None:
        """
        Creates a signature of the given fragment
        :param simulation: simulation containing the fragment
        :param fragment: molecule to be described by the signature generator
        :return: generated signature or None if the generation fails
        """
        if not self._is_allowed:
            return None

        try:
            return self._generate_signature(simulation, fragment)
        except SignatureGenerationError as exc:
            logger.info(f"[Signature] Signature generation of '{self.generator_name}' failed internally! {exc}")
        except Exception as exc:
            logger.error(f"[Signature] Signature generation of '{self.generator_name}' failed! {exc}")
            logger.exception(exc)

    def _generate_signature(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> str:
        """
        Implementation of generation of a signature of the given fragment
        :param simulation: simulation containing the fragment
        :param fragment: molecule to be described by the signature generator
        :raises SignatureGenerationError: if the generation fails
        :return: generated signature (on failure raises SignatureGenerationError)
        """
        raise NotImplementedError