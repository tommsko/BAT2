import json
from enum import Enum

from chembl_webresource_client.settings import Settings

Settings.Instance().TIMEOUT = 10
Settings.Instance().CACHING = False
Settings.Instance().TOTAL_RETRIES = 10
from .identifier_resolver import IdentifierResolver, IdentifierResolutionError
from ... import core_config

class NameDatabase(Enum):
    CCD = core_config.get("names", "ccd_db")
    MAPS  = core_config.get("names", "maps_db")
    MANUAL   = core_config.get("names", "manual_db")


class NameResolver(IdentifierResolver):
    """
    NameResolver attempts to resolve molecule names using various databases
    (no identification needs to be done)
    """

    def __init__(self, db_source: NameDatabase) -> None:
        """
        Initializes NameResolver for given signatures
        :param db_source: which database to use for name resolution
        """

        self.db_source: NameDatabase = db_source
        match db_source:
            case NameDatabase.CCD:
                resolver_name: str = "name_ccd"
                identifier_type: str = "INCHIKEY"
            case NameDatabase.MAPS:
                resolver_name: str = "name_maps"
                identifier_type: str = "INCHIKEY"
            case NameDatabase.MANUAL:
                resolver_name: str = "name_manual"
                identifier_type: str = "INCHIKEY"
            case _:
                raise RuntimeError("Unknown name database")

        super().__init__(resolver_name, identifier_type,["name"])


    def _is_applicable(self, signature: str) -> bool:
        """
        Checks whether signature is a name
        :param signature: textual signature to be identified
        :return: True if identifier name is applicable for signature, False otherwise
        """
        return True  # presupposing that output of name signature generator is good.
                     # No further checks are needed


    def _resolve_identity(self, signature: str, working_directory: str) -> list[str]:
        """
        Implementation for resolving signature on identity level
        :param signature: textual signature to be resolved
        :param working_directory: unused
        :raises IdentifierResolutionError: if resolution fails
        :return: list of identifiers exactly matched to the signature
        """

        with open(self.db_source.value, 'r') as fd:
            _db: dict = json.load(fd)

        if signature in _db:
            return [_db[signature]]

        return []


    def _resolve_similarity(self, signature: str, working_directory: str) -> list[tuple[str, float]]:
        """
        Implementation for resolving signature on similarity level
        :param signature: textual signature to be resolved
        :param working_directory: unused
        :raises IdentifierResolutionError: if resolution fails (no results are okay)
        :return: list of identifiers with associated similarity scores (0-1)
        """
        return [] # we are not doing name similarity resolution, it's a very bad idea


class NameResolverCCD(NameResolver):
    def __init__(self):
        """
        Initializes NAME identifier resolver w/ CCD name database (wrapper for NameResolver)
        """
        super().__init__(NameDatabase.CCD)

class NameResolverMAPS(NameResolver):
    def __init__(self):
        """
        Initializes NAME identifier resolver w/ MAPS name database (wrapper for NameResolver)
        """
        super().__init__(NameDatabase.MAPS)

class NameResolverMANUAL(NameResolver):
    def __init__(self):
        """
        Initializes NAME identifier resolver w/ MANUAL name database (wrapper for NameResolver)
        """
        super().__init__(NameDatabase.MANUAL)

