
from chembl_webresource_client.settings import Settings

Settings.Instance().TIMEOUT = 10
Settings.Instance().CACHING = False
Settings.Instance().TOTAL_RETRIES = 10
from chembl_webresource_client.new_client import new_client
from .identifier_resolver import IdentifierResolver, IdentifierResolutionError
from ... import core_config


class AtomResolver(IdentifierResolver):
    """
    AtomResolver just wraps identification process around elements identified by atom signature generator
    (no identification needs to be done)
    """

    def __init__(self) -> None:
        """
        Initializes AtomResolver for given signatures
        """
        super().__init__("atom", "ELEMENT",["element"])


    def _is_applicable(self, signature: str) -> bool:
        """
        Checks whether signature is atomic element
        :param signature: textual signature to be identified
        :return: True if identifier atom is applicable for signature, False otherwise
        """
        return True  # TODO: presupposing that output of atom signature generator is good.
                     # No further checks are needed


    def _resolve_identity(self, signature: str, working_directory: str) -> list[str]:
        """
        Implementation for resolving signature on identity level
        :param signature: textual signature to be resolved
        :param working_directory: unused
        :raises IdentifierResolutionError: if resolution fails
        :return: list of identifiers exactly matched to the signature
        """
        return [signature]


    def _resolve_similarity(self, signature: str, working_directory: str) -> list[tuple[str, float]]:
        """
        Implementation for resolving signature on similarity level
        :param signature: textual signature to be resolved
        :param working_directory: unused
        :raises IdentifierResolutionError: if resolution fails (no results are okay)
        :return: list of identifiers with associated similarity scores (0-1)
        """
        return [] # TODO: Atomic similarity is a smell, poignant one. But maybe

