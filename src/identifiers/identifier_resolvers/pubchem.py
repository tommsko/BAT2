import configparser
import time

import pubchempy
import requests

from .identifier_resolver import IdentifierResolver, IdentifierResolutionError
from ... import core_config


class PUBCHEMResolver(IdentifierResolver):
    """
    PUBCHEMResolver tries to resolve smiles identifiers using 2D structural PUBCHEM database search
    """

    def __init__(self) -> None:
        """
        Initializes PUBCHEMResolver for given signatures
        """
        super().__init__("pubchem_struct", "PUBCHEM_CID",["smiles"])


    def _is_applicable(self, signature: str) -> bool:
        """
        Checks whether signature is SMILES and applicable for selected PUBCHEM database
        :param signature: textual signature to be identified
        :return: True if identifier PUBCHEM is applicable for signature, False otherwise
        """
        return True  # presupposing that output of SMILES signature generator is good.
                     # No further checks are needed


    def _resolve_identity(self, signature: str, working_directory: str) -> list[str]:
        """
        Implementation for resolving signature on identity level
        :param signature: textual signature to be resolved
        :param working_directory: unused
        :raises IdentifierResolutionError: if resolution fails
        :return: list of identifiers exactly matched to the signature
        """
        return pubchem_fetch_identity_match(signature)


    def _resolve_similarity(self, signature: str, working_directory: str) -> list[tuple[str, float]]:
        """
        Implementation for resolving signature on similarity level
        :param signature: textual signature to be resolved
        :param working_directory: unused
        :raises IdentifierResolutionError: if resolution fails (no results are okay)
        :return: list of identifiers with associated similarity scores (0-1)
        """
        return pubchem_fetch_similarity_match(signature,
                                              threshold=core_config.getint("pubchem",
                                                                           "similarity_threshold")
            )


def _pubchem_cids_to_inchikeys(cids: list[int]) -> list[str]:
    """
    Translates PUBCHEM's molecular CIDs to InChI Keys using PUBCHEM servers (using API)
    :param cids: list of cids to be translated
    :return: equivalent list of InChIs
    """
    inchi_keys: list[str] = []
    for cid in cids:
        if cid == 0:
            continue  # CID is never zero, just a result of bad POST fallback request
        try:
            inchi_keys.append(pubchempy.get_properties("InChIKey", cid)[0]["InChIKey"])
        except Exception:
            ...  # resolution can fail, their databases are not perfect
    return inchi_keys


def _pubchem_fetch_identity_match_post(smiles: str) -> list[int]:
    """
    As a fallback, attempts direct (POST) search of PUBCHEM's servers for identity match of SMILES
    :param smiles: signature of molecule to search
    :raises IdentifierResolutionError on any failure
    :return: list of CID's
    """
    try:
        cids_text: str = requests.post(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/txt",
            data={"smiles": smiles},
        ).text.strip()
        if "status" in cids_text.lower():  # request failed
            return []
        return [int(cid) for cid in cids_text.split("\n") if cid.isnumeric()]
    except Exception as exc:
        if 'ServerBusy' in str(exc):
            return []  # we don't want to fail the whole search just because server is busy
        raise IdentifierResolutionError(
            f"[Identifier] PUBCHEM POST identity search failed: {exc}"
        )


def pubchem_fetch_identity_match(smiles: str) -> list[int]:
    """
    Searches PUBCHEM's databases for identity match of SMILES
    :param smiles: identifier of molecule to search
    :raises IdentifierResolutionError: on PUBCHEM search failure (not finding anything is NOT a failure)
    :return: list of InChI Keys of identical molecules
    """

    retry_attempts: int = core_config.getint("pubchem", "retry_attempts")
    retry_delay: int = core_config.getint("pubchem", "retry_delay_s")

    for _ in range(retry_attempts):
        try:
            try:
                exact_matches: list[pubchempy.Compound] = pubchempy.get_compounds(
                    smiles, namespace="smiles", searchtype="identity", as_dataframe=False
                )
                return [int(match.cid) for match in exact_matches]
            except pubchempy.BadRequestError:  # fallback
                return _pubchem_fetch_identity_match_post(smiles)
        except pubchempy.TimeoutError:  # this is expected
            time.sleep(retry_delay)
        except Exception as exc:
            if 'ServerBusy' in str(exc):  # this is also expected
                time.sleep(retry_delay)
            else:
                raise IdentifierResolutionError(
                    f"[Identifier] PUBCHEM API identity search failed: {exc}"
                )
    return []


def _pubchem_fetch_similar_match_post(smiles: str, threshold: int) -> list[int]:
    """
    As a fallback, attempts direct (POST) search of PUBCHEM's servers for similarity match of SMILES
    :param smiles: identifier of molecule to search
    :param threshold: threshold of similarity
    :raises IdentifierResolutionError on any failure
    :return: list of CID's
    """
    try:
        cids_text: str = requests.post(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/cids/txt?Threshold={threshold}",
            data={"smiles": smiles},
        ).text.strip()
        if "status" in cids_text.lower():  # request failed
            return []
        return [int(cid) for cid in cids_text.split("\n") if cid.isnumeric()]
    except Exception as exc:
        if 'ServerBusy' in str(exc):
            return []  # we don't want to fail the whole search just because server is busy
        raise IdentifierResolutionError(
            f"[Identifier] PUBCHEM POST similarity search failed: {exc}"
        )


def pubchem_fetch_similarity_match(smiles: str, threshold: int) -> list[tuple[int, float]]:
    """
    Searches PUBCHEM's databases for similar match of SMILES identifier given the threshold.
    Low threshold will return more matches, but also it may result in timeout given PUBCHEM servers
    :param smiles: identifier of molecule to search
    :param threshold: minimum similarity threshold for filtering hits on PUBCHEM side
    :raises IdentifierResolutionError: on PUBCHEM search failure (excluding timeout)
    :return: list of InChI Keys of identical molecules and similarity percentage they were found with
    """

    retry_attempts: int = core_config.getint("pubchem", "retry_attempts")
    retry_delay: int = core_config.getint("pubchem", "retry_delay_s")

    for _ in range(retry_attempts):
        try:
            try:
                compounds: list[pubchempy.Compound] = pubchempy.get_compounds(
                    smiles,
                    namespace="smiles",
                    searchtype="similarity",
                    as_dataframe=False,
                    Threshold=threshold,
                )
                cids: list[int] = [int(match.cid) for match in compounds]
            except pubchempy.BadRequestError:  # fallback
                cids: list[int] = _pubchem_fetch_similar_match_post(smiles, threshold)
            return list(
                zip(
                    cids,
                    [threshold for _ in range(len(cids))],
                )
            )
        except pubchempy.TimeoutError:  # this is expected
            time.sleep(retry_delay)
        except Exception as exc:
            if 'ServerBusy' in str(exc):  # this is also expected
                time.sleep(retry_delay)
            else:
                raise IdentifierResolutionError(
                    f"[Identifier] PUBCHEM API similarity search failed: {exc}"
                )
    return []
