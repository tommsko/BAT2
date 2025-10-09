import configparser
import json
import logging
import os.path
import subprocess
import time
import requests
from requests import Response
from ... import core_config
from .identifier_resolver import IdentifierResolver, IdentifierResolutionError
from ...utils.file_utils import mk_random_filename
log = logging.getLogger("identifier")

class AlphaFindStructResolver(IdentifierResolver):
    """
    AlphaFindStructResolver tries to protein/nucleic acid structure using AlphaFind structure database search
    """

    def __init__(self) -> None:
        """
        Initializes AlphaFindStructResolver for given signatures
        """
        super().__init__("alphafind_struct", "ALPHAFOLD_ID",["pdb_protein"])


    def _is_applicable(self, signature: str) -> bool:
        """
        Checks whether signature is PDB-like and applicable for selected AlphaFind database
        :param signature: textual signature to be identified
        :return: True if identifier AlphaFindStructResolver is applicable for signature, False otherwise
        """
        return True  # TODO: presupposing that output of pdb signature generator is good. No further checks are needed


    def _resolve_identity(self, signature: str, working_directory: str) -> list[str]:
        """
        Implementation for resolving signature on identity level
        :param signature: textual signature to be resolved
        :param working_directory: working directory of the analysis
        :raises IdentifierResolutionError: if resolution fails
        :return: list of identifiers exactly matched to the signature
        """

        if core_config.getboolean("miscellaneous", "use_temporary_directory"):
            tmp_dir: str = core_config.get("miscellaneous", "temporary_directory")
        else:
            tmp_dir: str = working_directory

        pdb_path: str = mk_random_filename("pdb",
                                           tmp_dir)
        with open(pdb_path, 'w') as fd:
            fd.write(signature)

        try:
            return alphafind_identity(pdb_path)
        except Exception as exc:
            log.warning(f"[Identifier] PDB structure resolution failed: {exc}")
            return []
        finally:
            os.unlink(pdb_path)
            # will not overwrite return
            # https://stackoverflow.com/questions/19805654/python-try-finally-block-returns



    def _resolve_similarity(self, signature: str, working_directory: str) -> list[tuple[str, float]]:
        """
        Implementation for resolving signature on similarity level
        :param signature: textual signature to be resolved
        :param working_directory: working directory of the analysis
        :raises IdentifierResolutionError: if resolution fails (no results are okay)
        :return: list of identifiers with associated similarity scores (0-1)
        """

        if core_config.getboolean("miscellaneous", "use_temporary_directory"):
            tmp_dir: str = core_config.get("miscellaneous", "temporary_directory")
        else:
            tmp_dir: str = working_directory


        pdb_path: str = mk_random_filename("pdb",
                                           tmp_dir)
        with open(pdb_path, 'w') as fd:
            fd.write(signature)

        try:
            return alphafind_similarity(pdb_path)
        except Exception as exc:
            log.warning(f"[Identifier] PDB structure resolution failed: {exc}")
            return []
        finally:
            os.unlink(pdb_path)
            # will not overwrite return
            # https://stackoverflow.com/questions/19805654/python-try-finally-block-returns


def _alphafind_curl(pdb_path: str, request_timeout: int) -> str:
    """
    Executes AlphaFind CURL request
    :param pdb_path: to pdb file
    :param request_timeout: maximum number of seconds for request
    :return: string response
    """
    try:
        return subprocess.check_output(
            ['curl',
             '--request', 'POST',
             '--url', 'https://api.stage.alphafind.dyn.cloud.e-infra.cz/search?limit=10',
             #'--http1.1',
             '--header', 'content-type: multipart/form-data',
             '--form', f'file=@{pdb_path}'], timeout=request_timeout).decode()
    except Exception:
        return ""


def _alphafind_fetch_response(pdb_path: str) -> dict | None:
    """
    Executes remote AlphaFind similarity search
    :param pdb_path: of local .pdb file to be searched
    :return: json response if successful, None on timeout
    """
    log.debug(f"[Identifier] AlphaFind resolving similarity search for '{pdb_path}'")

    retry_count: int = core_config.getint("ALPHAFIND", 'retry_count')
    retry_delay: int = core_config.getint("ALPHAFIND", 'retry_delay_s')
    request_timeout: int = core_config.getint("ALPHAFIND", 'request_timeout_s')

    _alphafind_curl(pdb_path, request_timeout)  # initial request
    log.debug("AlphaFind CURL request initiated")

    for _ in range(retry_count):
        time.sleep(retry_delay)
        try:
            log.debug("AlphaFind CURL request check")
            response: dict = json.loads(_alphafind_curl(pdb_path, request_timeout))
            log.debug(f"... {response}")
            if 'results' in response:
                return response
        except Exception:
            ...

    return None


def alphafind_similarity(pdb_path: str) -> list[tuple[str, float]]:
    """
    Executes AlphaFind similarity search
    :param pdb_path: of local .pdb file to be searched
    :return: list of alphafold identifiers with associated aligned_percentage
    """

    response: dict | None = _alphafind_fetch_response(pdb_path)
    if response is None or 'results' not in response:
        return []

    results: list[tuple[str, float]] = []
    for entry in response['results']:
        if 'object_id' in entry and 'aligned_percentage' in entry:
            results.append((entry['object_id'], entry['aligned_percentage']))
    return results


def alphafind_identity(pdb_path: str) -> list[str]:
    """
    Executes AlphaFind identity search
    :param pdb_path: of local .pdb file to be searched
    :return: list of alphafold identifiers
    """

    return [entry[0] for entry in alphafind_similarity(pdb_path) if entry[1] == 1.0]
