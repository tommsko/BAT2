import configparser
import json
import logging
import os.path
import subprocess
import requests
from requests import Response

from ... import core_config
from .identifier_resolver import IdentifierResolver, IdentifierResolutionError
from ...utils.file_utils import mk_random_filename

log = logging.getLogger("identifier")

class PDBStructResolver(IdentifierResolver):
    """
    PDBStructResolver tries to protein/nucleic acid structure using PDB structure database search
    """

    def __init__(self, as_protein: bool, as_nucleotide: bool) -> None:
        """
        Initializes PDBStructResolver for given signatures
        """
        assert as_protein or as_nucleotide, "One mode of operation is necessary"
        self._protein_mode: bool = as_protein
        self._nucleotide_mode: bool = as_nucleotide

        super().__init__("pdb_struct_protein" if as_protein else "pdb_struct_nucleotide",
                         "PDB_ID",
                         ["pdb_protein" if as_protein else "pdb_nucleic"])


    def _is_applicable(self, signature: str) -> bool:
        """
        Checks whether signature is PDB-like and applicable for selected PDB database
        :param signature: textual signature to be identified
        :return: True if identifier PDBStructResolver is applicable for signature, False otherwise
        """
        return True  # presupposing that output of pdb signature generator is good. No further checks are needed


    def _resolve_identity(self, signature: str, working_directory: str) -> list[str]:
        """
        Implementation for resolving signature on identity level
        :param signature: textual signature to be resolved
        :param working_directory: working directory of the analysis
        :raises IdentifierResolutionError: on resolution failure
        :return: list of identifiers exactly matched to the signature
        """
        if core_config.getboolean("miscellaneous", "use_temporary_directory"):
            tmp_dir: str = core_config.get("miscellaneous", "temporary_directory")
        else:
            tmp_dir: str = working_directory

        pdb_path: str = mk_random_filename("pdb", tmp_dir)
        with open(pdb_path, 'w') as fd:
            fd.write(signature)

        try:
            deployment_file_access: bool = core_config.getboolean("PDB", "use_deployment_file_access")
            return pdb_fetch_identical_structures(pdb_path, deployment_file_access, self._protein_mode, self._nucleotide_mode)
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
            deployment_file_access: bool = core_config.getboolean("PDB", "use_deployment_file_access")
            return pdb_fetch_similar_structures(pdb_path, deployment_file_access, self._protein_mode, self._nucleotide_mode)
        except Exception as exc:
            log.warning(f"[Identifier] PDB structure resolution failed: {exc}")
            return []
        finally:
            os.unlink(pdb_path)
            # will not overwrite return
            # https://stackoverflow.com/questions/19805654/python-try-finally-block-returns

class PDBStructResolverProtein(PDBStructResolver):
    def __init__(self):
        super().__init__(as_protein=True, as_nucleotide=False)

class PDBStructResolverNucleotide(PDBStructResolver):
    def __init__(self):
        super().__init__(as_protein=False, as_nucleotide=True)

def _replace_placeholders(command: list[str], placeholders: dict[str, str]) -> None:
    """
    Replaces placeholders with values in command to be executed
    :param command: to replace placeholders in (structured as input of subprocess.run)
    :param placeholders: mapping of string replacements to be done
    :return: None (modifies command parameter)
    :raises KeyError: if not all placeholders were replaced
    """
    placeholders_to_replace: set[str] = set(placeholders.keys())
    for i in range(len(command)):
        command_part: str = command[i]
        for placeholder_from, placeholder_to in placeholders.items():
            if placeholder_from in command_part:
                command[i] = command[i].replace(placeholder_from, placeholder_to)
                placeholders_to_replace.remove(placeholder_from)

    if len(placeholders_to_replace) != 0:
        raise KeyError(
            f"No substitution possible for placeholders: {placeholders_to_replace}"
        )


def _upload_file(pdb_path: str) -> str | None:
    """
    Uploads PDB file to remote server
    :param pdb_path: to the .pdb file
    :return: name of the .pdb file, if successful
    :raises CalledProcessError: on upload error
    """
    pdb_name: str = os.path.basename(pdb_path)
    log.debug(f"[Identifier] PDB uploading '{pdb_name}'")

    upload_command: list[str] = [
        _v for _, _v in core_config.items("PDB_UPLOAD_UPLOAD_COMMAND")
    ]
    _replace_placeholders(
        upload_command, {"{FILEPATH}": pdb_path, "{FILENAME}": pdb_name}
    )
    subprocess.run(upload_command, check=True)
    log.debug("[Identifier] PDB upload successful")
    return pdb_name


def _delete_file(pdb_name: str) -> None:
    """
    Removes PDB file from remote server
    :param pdb_name: path to the .pdb file
    :return: None
    :raises Exception: on ssh error
    """
    log.debug(f"[Identifier] PDB deleting '{pdb_name}'")
    purge_command: list[str] = [
        _v for _, _v in core_config.items("PDB_UPLOAD_DELETE_COMMAND")
    ]
    _replace_placeholders(purge_command, {"{FILENAME}": pdb_name})
    subprocess.run(purge_command, check=True)
    log.debug("[Identifier] PDB delete successful")


def _fetch_pdb_similarity(pdb_name: str, deployment_file_access: bool, as_protein: bool, as_nucleic: bool) -> dict:
    """
    Executes online PDB structural similarity search given uploaded file
    :param pdb_name: name of the PDB file
    :param deployment_file_access: if set to True, file should be accessible using API
    :param as_protein: search protein entries
    :param as_nucleic: search nucleic acid entries
    :return: json dictionary of search results
    :raises Exception: on POST request error
    """
    if deployment_file_access:
        STRUCT_URL = core_config.get("PDB", "deployment_file_url")
        assert "{FILENAME}" in STRUCT_URL, "illegal template for URL"
        assert "{UUID}" in STRUCT_URL, "illegal template for URL"
        remote_path: str = STRUCT_URL.replace("{FILENAME}", pdb_name)
        remote_path = remote_path.replace("{UUID}", core_config.get("TMP", "CURR_UUID"))
    else:
        STRUCT_URL = core_config.get("PDB_UPLOAD_URL", "url")
        assert "{FILENAME}" in STRUCT_URL, "illegal template for URL"
        remote_path: str = STRUCT_URL.replace("{FILENAME}", pdb_name)

    url = "https://search.rcsb.org/rcsbsearch/v2/query"

    required_entity_type: str = "Protein (only)" if as_protein else "Nucleic acid (only)"
    request = {
          "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
              {
                "type": "terminal",
                "service": "structure",
                "parameters": {
                  "value": { "url": remote_path,"format": "pdb"},
                  "operator": "relaxed_shape_match"
                }
              },
              {
                "type": "terminal",
                "service": "text",
                "parameters": {
                  "attribute": "rcsb_entry_info.selected_polymer_entity_types",
                  "operator": "in",
                  "value": [required_entity_type]
                }
              }
            ]
          },
          "return_type": "entry"
        }

    log.debug(f"[Identifier] [DB querying similarity for '{remote_path}'")
    response: Response = requests.post(url, json=request)
    response_txt: str = response.text
    response.close()
    response_json: dict = json.loads(response_txt)
    return response_json


def _parse_pdb_response(response: dict) -> list[tuple[str, float]]:
    """
    Parses PDB similarity search response and returns list of identifiers and their similarity percentages
    :param response: result of _fetch_pdb_similarity(...)
    :return: list of identifiers and their similarity percentages
    """
    if "result_set" not in response:
        return []

    results: list[tuple[str, float]] = []
    for hit in response["result_set"]:
        results.append((hit["identifier"], hit["score"] * 100))
    return results


def pdb_fetch_similar_structures(pdb_path: str, deployment_file_access: bool, as_protein: bool, as_nucleic: bool) -> list[tuple[str, float]]:
    """
    Executes remote PDB similarity search
    :param pdb_path: of local .pdb file to be searched
    :param deployment_file_access: if set to True, file should be available through API without uploading
    :param as_protein: search protein entries
    :param as_nucleic: search nucleic acid entries
    :return: list of PDB identifiers as well as similarity percentages
    """
    assert as_protein or as_nucleic, "One mode of operation is necessary"

    log.debug(f"[Identifier] PDB resolving similarity search '{pdb_path}'")

    if deployment_file_access:
        filename: str = os.path.basename(pdb_path)
    else:
        filename: str = _upload_file(pdb_path)

    response: dict = _fetch_pdb_similarity(filename, deployment_file_access, as_protein=as_protein, as_nucleic=as_nucleic)

    similar_proteins: list[tuple[str, float]] = _parse_pdb_response(response)
    if not deployment_file_access:
        _delete_file(filename)
    return similar_proteins


def pdb_fetch_identical_structures(pdb_path: str, deployment_file_access: bool, as_protein: bool, as_nucleic: bool) -> list[str]:
    """
    Executes remote PDB identity search
    :param pdb_path: of local .pdb file to be searched
    :param deployment_file_access: if set to True, file should be available through API without uploading
    :param as_protein: search protein entries
    :param as_nucleic: search nucleic acid entries
    :return: list of PDB identifiers
    """

    assert as_protein or as_nucleic, "One mode of operation is necessary"

    log.debug(f"resolving PDB identity search '{pdb_path}' (protein={as_protein}, nucleic={as_nucleic})")
    return [
        pdb_id
        for pdb_id, similarity in pdb_fetch_similar_structures(pdb_path, deployment_file_access, as_protein=as_protein, as_nucleic=as_nucleic)
        if similarity == 100
    ]
