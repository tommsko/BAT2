import configparser
import logging
import os.path
import sqlite3
import warnings
from .. import core_config


log = logging.getLogger("base")

def create_empty_cache_if_not_exists(cache_path: str) -> None:
    """
    If cache does not exist, creates it and populates it with required tables
    :param cache_path: to be used
    :return: None
    """
    if os.path.exists(cache_path):
        return

    log.warning(f"[Cache] Resolver cache not found at '{cache_path}': Initializing empty")

    directory_path: str = os.path.dirname(cache_path)
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

    with sqlite3.connect(cache_path) as conn:
        conn.execute(
            """CREATE TABLE "CACHE_IDENTITY" (
                        "RESOLVER_NAME"	TEXT NOT NULL,
                        "SIGNATURE"	    TEXT NOT NULL,
                        "IDENTIFIER"	TEXT NOT NULL,
                        PRIMARY KEY("IDENTIFIER","SIGNATURE","RESOLVER_NAME"));
                    """
        )

        conn.execute(
            """CREATE TABLE "CACHE_SIMILARITY" (
                        "RESOLVER_NAME"	TEXT NOT NULL,
                        "SIGNATURE"	    TEXT NOT NULL,
                        "IDENTIFIER"	TEXT NOT NULL,
                        "SIMILARITY"	REAL,
                        PRIMARY KEY("IDENTIFIER","SIGNATURE","RESOLVER_NAME"));
                    """
        )

        conn.commit()


class ResolverCache:
    def __init__(self) -> None:
        """
        Initializes resolver cache client (and creates database if it does not exist)
        """

        cache_directory: str = core_config.get("identifiers_cache", "directory")
        cache_file: str      = core_config.get("identifiers_cache", "file")
        cache_path: str      = os.path.join(cache_directory, cache_file)

        create_empty_cache_if_not_exists(cache_path)
        self.conn: sqlite3.Connection = sqlite3.connect(cache_path)

    def save_identity_identifiers(self, resolver_name: str, signature: str, identifiers: list[str]) -> None:
        """
        Saves all (new) identifier associated with the signature by exact-match into the cache
        Old identifiers are not overwritten/deleted
        :param resolver_name: of the resolver that found the identifiers
        :param signature: for which identifiers were found
        :param identifiers: to be saved
        :return: None
        """

        with self.conn as conn:
            identifiers_already_cached: set[str] = set(
                self.fetch_identity_identifiers(resolver_name, signature)
            )
            for ident in identifiers:
                if ident in identifiers_already_cached:
                    continue
                conn.execute(
                    "INSERT INTO CACHE_IDENTITY VALUES (?, ?, ?)",
                    (resolver_name, signature, ident),
                )
            conn.commit()

    def fetch_identity_identifiers(
        self, resolver_name: str, signature: str
    ) -> list[str]:
        """
        Extracts all exact-match identifiers associated with given signature and resolver from cache
        :param resolver_name: of resolver trying to identify the signature
        :param signature: for which to fetch exact-match identifiers
        :return: list of all identifiers associated with the resolver and signature on exact-match basis
        """
        with self.conn as conn:
            records = conn.execute(
                "SELECT IDENTIFIER FROM CACHE_IDENTITY WHERE RESOLVER_NAME=? AND SIGNATURE=?",
                (resolver_name, signature),
            ).fetchall()
            return [record[0] for record in records]

    def save_similarity_identifiers(
        self, resolver_name: str, signature: str, identifiers: list[tuple[str, float]]
    ) -> None:
        """
        Saves all (new) identifier associated with the signature by similarity search into the cache
        Old identifiers are not overwritten/deleted
        :param resolver_name: of the resolver that found the identifiers
        :param signature: for which identifiers were found
        :param identifiers: to be saved
        :return: None
        """

        with self.conn as conn:
            identifiers_already_cached: set[str] = set(
                _ident
                for _ident, _ in self.fetch_similarity_identifiers(
                    resolver_name, signature
                )
            )
            for _ident, _similarity in identifiers:
                if _ident in identifiers_already_cached:
                    continue
                conn.execute(
                    "INSERT INTO CACHE_SIMILARITY VALUES (?, ?, ?, ?)",
                    (resolver_name, signature, _ident, _similarity),
                )
            conn.commit()

    def fetch_similarity_identifiers(
        self, resolver_name: str, signature: str
    ) -> list[tuple[str, float]]:
        """
        Extracts all similarity-based identifiers associated with given signature and resolver from cache
        :param resolver_name: of resolver trying to identify the signature
        :param signature: for which to fetch similarity-based identifiers
        :return: list of all identifiers associated with the resolver and signature on similarity basis
        """
        with self.conn as conn:
            records = conn.execute(
                "SELECT IDENTIFIER, SIMILARITY FROM CACHE_SIMILARITY "
                "WHERE RESOLVER_NAME=? AND SIGNATURE=?",
                (resolver_name, signature),
            ).fetchall()
            return records




