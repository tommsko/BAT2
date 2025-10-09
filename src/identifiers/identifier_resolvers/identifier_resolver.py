import logging
from logging import Logger

from ... import core_config
from ...loader.simulation import Simulation, SegmentSupports
from ..cache_manager import ResolverCache

logger: Logger = logging.getLogger("identifier")

class IdentifierResolutionError(Exception):
    def __init__(self, message):
        super().__init__(message)


class IdentifierResolver:
    """
    Abstract class for identifier resolvers which input signature and output a specific identifier (database id string)
    """
    def __init__(self, resolver_name: str, identifier_type: str, processes_signature_types: list[str]) -> None:
        """
        Initializes IdentifierResolver for given signatures
        :param resolver_name: name of the resolver
        :param identifier_type: string type of identifier (if resolved)
        :param processes_signature_types: list of signature types (signature_generator.py::__init__::signature_type)
                                          processed by the resolver
        """
        self.resolver_name: str = resolver_name
        self.identifier_type: str = identifier_type
        self.accepted_signatures: list[str] = processes_signature_types

        self.allowed_identity: bool = core_config.getboolean("identifiers",
                                                             self.resolver_name,
                                                             fallback=False)
        self.allowed_similarity: bool = core_config.getboolean("identifiers",
                                                               self.resolver_name + '_similarity',
                                                               fallback=False)

    def is_applicable(self, signature: str, signature_type: str, as_similarity: bool) -> bool:
        """
        Checks whether identifier resolver is applicable for given signature
        :param signature: textual signature to be identified
        :param signature_type: type of signature to be identified
        :param as_similarity: if True, checks for similarity resolution instead of identity
        :return: True if identifier resolver is applicable for signature and its type
        """
        if not self.allowed_identity or (as_similarity and not self.allowed_similarity):
            return False

        if signature_type not in self.accepted_signatures:
            return False

        return self._is_applicable(signature)


    def _is_applicable(self, signature: str) -> bool:
        """
        Implementation of check whether the identifier resolver is applicable to the given signature
        :param signature: textual signature to be identified
        :return: True if identifier resolver is applicable for signature, False otherwise
        """
        raise NotImplementedError


    def resolve_signature(self, signature: str, as_similarity: bool, working_directory: str) -> list[str] | list[tuple[str, float]] | None:
        """
        Attempts to resolve signature to identifier (database entry ID)
        :param signature: to be resolved
        :param as_similarity: if True, similarity search is used instead of exact-match
        :param working_directory: working directory of the analysis
        :return: list of resolved identifier in case of exact-match
                 or list of identifiers with associated similarity scores (0-1) in case of similarity
                 or None if identification fails / is not allowed
        """
        if (not as_similarity and not self.allowed_identity) or (as_similarity and not self.allowed_similarity):
            return None

        cache: ResolverCache = ResolverCache()

        try:
            if not as_similarity:
                if (cached := cache.fetch_identity_identifiers(self.resolver_name, signature)):
                    return cached

                idents: list[str] = self._resolve_identity(signature, working_directory)
                if len(idents) == 0:
                    return None
                cache.save_identity_identifiers(self.resolver_name, signature, idents)
                return idents

            else:
                if (cached := cache.fetch_similarity_identifiers(self.resolver_name, signature)):
                    return cached

                idents: list[tuple[str, float]] = self._resolve_similarity(signature, working_directory)
                if len(idents) == 0:
                    return None
                cache.save_similarity_identifiers(self.resolver_name, signature, idents)
                return idents
        except IdentifierResolutionError as exc:
            logger.info(f"[Identifier] Identifier resolution of '{self.resolver_name}' failed internally! {exc}")
        except Exception as exc:
            logger.error(f"[Identifier] Identifier resolution of '{self.resolver_name}' failed! {exc}")
            logger.exception(exc)

        return None


    def _resolve_identity(self, signature: str, working_directory: str) -> list[str]:
        """
        Implementation for resolving signature on identity level
        :param signature: textual signature to be resolved
        :param working_directory: working directory of the analysis
        :raises IdentifierResolutionError: if resolution fails or nothing is found
        :return: list of identifiers exactly matched to the signature
        """
        raise NotImplementedError


    def _resolve_similarity(self, signature: str, working_directory: str) -> list[tuple[str, float]]:
        """
        Implementation for resolving signature on similarity level
        :param signature: textual signature to be resolved
        :param working_directory: working directory of the analysis
        :raises IdentifierResolutionError: if resolution fails (no results are okay)
        :return: list of identifiers with associated similarity scores (0-1)
        """
        raise NotImplementedError

