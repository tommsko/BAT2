import json
import os
import subprocess
from dataclasses import dataclass
from enum import Enum
from typing import Callable

from ... import core_config
from ...constants import KNOWN_RESIDUE_CODES, KNOWN_NUCLEIC_ACIDS_CODES
from .identifier_resolver import IdentifierResolver, IdentifierResolutionError
from ...utils.file_utils import mk_random_filename

class BlastDatabase(Enum):
    PDB_PROT = core_config.get("blast", "pdb_aa_db")
    UNIPROT  = core_config.get("blast", "uniprot_aa_db")
    PDB_NT   = core_config.get("blast", "pdb_nt_db")


class BLASTResolver(IdentifierResolver):
    """
    BLASTResolver tries to identify residue sequences (amino acids or nucleotides) using various databases
    """

    def __init__(self, db_source: BlastDatabase) -> None:
        """
        Initializes BLASTResolver for given signatures
        :param db_source: which database to use for blast
        """

        is_protein: bool = db_source in (BlastDatabase.PDB_PROT, BlastDatabase.UNIPROT)

        self.db_source: BlastDatabase = db_source
        match db_source:
            case BlastDatabase.PDB_PROT:
                resolver_name: str = "blast_pdb_protein"
                identifier_type: str = "PDB_ID"
            case BlastDatabase.UNIPROT:
                resolver_name: str = "blast_uniprot"
                identifier_type: str = "UNIPROT_ID"
            case BlastDatabase.PDB_NT:
                resolver_name: str = "blast_pdb_nucleic"
                identifier_type: str = "PDB_ID"
            case _:
                raise RuntimeError("Unknown blast database")

        super().__init__(resolver_name, identifier_type,
                         ["residue_sequence_protein" if is_protein
                                                else "residue_sequence_nucleotide"])


    def _is_applicable(self, signature: str) -> bool:
        """
        Checks whether signature is fasta-like (fasta without header) and applicable for selected blast database
        :param signature: textual signature to be identified
        :return: True if identifier BLASTResolver with selected database is applicable for signature, False otherwise
        """
        if self.db_source in (BlastDatabase.PDB_PROT, BlastDatabase.UNIPROT):
            return all(c in KNOWN_RESIDUE_CODES.values() for c in signature)

        return all(c in KNOWN_NUCLEIC_ACIDS_CODES.values() for c in signature)

    def _get_exec_path(self) -> str:
        """
        Retrieves BLAST executable path from configuration given the selected database
        :return: path to BLAST executable
        """
        if self.db_source in (BlastDatabase.PDB_PROT, BlastDatabase.UNIPROT):
            return core_config.get("blast", "blastp_exec")
        else:
            return core_config.get("blast", "blastn_exec")

    def _resolve_identity(self, signature: str, working_directory: str) -> list[str]:
        """
        Implementation for resolving signature on identity level
        :param signature: textual signature to be resolved
        :param working_directory: working directory of the analysis
        :raises IdentifierResolutionError: if resolution fails
        :return: list of identifiers exactly matched to the signature
        """
        return run_blast_identity(signature,
                                  self._get_exec_path(),
                                  self.db_source.value,
                                  self.db_source == BlastDatabase.PDB_NT,
                                  working_directory)

    def _resolve_similarity(self, signature: str, working_directory: str) -> list[tuple[str, float]]:
        """
        Implementation for resolving signature on similarity level
        :param signature: textual signature to be resolved
        :param working_directory: working directory of the analysis
        :raises IdentifierResolutionError: if resolution fails (no results are okay)
        :return: list of identifiers with associated similarity scores (0-1)
        """
        return run_blast_similarity(signature,
                                  self._get_exec_path(),
                                  self.db_source.value,
                                  self.db_source == BlastDatabase.PDB_NT,
                                  working_directory)

class BLASTResolverPDBAA(BLASTResolver):
    def __init__(self):
        """
        Initializes BLAST identifier resolver w/ PDB protein database (wrapper for BLASTResolver)
        """
        super().__init__(BlastDatabase.PDB_PROT)

class BLASTResolverPDBNT(BLASTResolver):
    def __init__(self):
        """
        Initializes BLAST identifier resolver w/ PDB nucleotide database (wrapper for BLASTResolver)
        """
        super().__init__(BlastDatabase.PDB_NT)

class BLASTResolverUNIPROT(BLASTResolver):
    def __init__(self):
        """
        Initializes BLAST identifier resolver w/ UNIPROT protein database (wrapper for BLASTResolver)
        """
        super().__init__(BlastDatabase.UNIPROT)

def _run_blast_local(sequence: str,
                     executable_path: str,
                     database_path: str,
                     is_nucleotide: bool,
                     working_directory: str) -> list[dict] | None:
    """
    Executes BLAST search with local engine
    :param sequence: to be searched
    :param executable_path: path to appropriate (protein or nucleotide) BLAST search engine
    :param database_path: path to BLAST database
    :param is_nucleotide: if True, utilizes nucleotide optimization parameters (instead of protein ones)
    :param working_directory: working directory of the analysis
    :return: list of jsons (results) returned, or None on failure
    :raises: IdentifierResolutionError on blast failure
    """

    if core_config.getboolean("miscellaneous", "use_temporary_directory"):
        tmp_dir: str = core_config.get("miscellaneous", "temporary_directory")
    else:
        tmp_dir: str = working_directory

    output_path: str = mk_random_filename("zip",
                                          tmp_dir)

    try:
        if len(sequence) < core_config.getint('blast',
                                              'local_instance_run_short_search_if_under_X_residues'):
            short_search: list[str] = ['-task', 'blastp-short'] if not is_nucleotide else ['-task', 'blastn-short']
        else:
            short_search: list[str] = []

        result = subprocess.run(
            [
                executable_path,
                "-outfmt",
                "15",
                "-db",
                database_path,
                "-out",
                output_path,
            ] + short_search,
            input=sequence,
            text=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        assert result.returncode == 0, "non-zero exit code"
    except Exception as exc:
        raise IdentifierResolutionError(f"Local BLAST failed: {exc}")

    with open(output_path, "r") as fd:
        json_files: list[dict] = [json.load(fd)]

    os.unlink(output_path)
    return json_files


@dataclass
class BlastHit:
    accessions: list[str]
    evalue: float
    bits: float
    identity_p: float
    coverage_p: float


def _process_blast_results(json_files: list[dict]) -> list[BlastHit]:
    """
    Parses the BLAST results and standardizes the format of hits
    :param json_files: output of _run_blast*
    :return: list of hits
    """

    if json_files is None:
        return []

    results: list[BlastHit] = []

    for data in json_files:
        if "BlastOutput2" not in data:
            continue

        blast: dict = data["BlastOutput2"]
        if isinstance(blast, list):
            blast = blast[0]

        if (
            "report" not in blast
            or "results" not in blast["report"]
            or "search" not in blast["report"]["results"]
            or "hits" not in blast["report"]["results"]["search"]
        ):
            continue

        for hit in blast["report"]["results"]["search"]["hits"]:
            accessions: list[str] = [match["accession"] for match in hit["description"]]

            query_len: int = blast["report"]["results"]["search"]["query_len"]
            #hit_len: int = hit["len"]

            _scores_section: dict = hit["hsps"][0]
            bits: float = _scores_section["bit_score"]
            evalue: float = _scores_section["evalue"]
            aligned_residues: int = _scores_section["align_len"]
            identity_p: float = _scores_section["identity"] / aligned_residues
            coverage_p: float = aligned_residues / query_len
            results.append(BlastHit(accessions, evalue, bits, identity_p, coverage_p))
    return results


def _filter_blast_results(
    results: list[BlastHit], by=Callable[[BlastHit], bool]
) -> list[BlastHit]:
    """
    Simple filter wrapper for BLAST results
    :param results: BLAST hits to filter
    :param by: lambda defining the filter (keep) condition
    :return: BLAST hits following the filter
    """
    return [hit for hit in results if by(hit)]


def run_blast_identity(sequence: str,
                       executable_path: str,
                       database_path: str,
                       is_nucleotide: bool,
                       working_directory: str) -> list[str]:
    """
    Executes identity search with BLAST
    :param sequence: to be searched
    :param executable_path: path to appropriate (protein or nucleotide) BLAST search engine
    :param database_path: path to BLAST database
    :param is_nucleotide: if True, utilizes nucleotide optimization parameters (instead of protein ones)
    :param working_directory: working directory of the analysis
    :return: accessions of identical proteins
    :raises IdentifierResolutionError: on BLAST failure (finding no results is NOT a failure)
    """

    hits: list[BlastHit] = _process_blast_results(_run_blast_local(sequence, executable_path,
                                                                   database_path, is_nucleotide, working_directory))

    hits = _filter_blast_results(
        hits,
        lambda hit: hit.evalue <= core_config.getfloat("blast", "exact_match_max_evalue"),
    )
    hits = _filter_blast_results(
        hits, lambda hit: hit.bits >= core_config.getint("blast", "exact_match_min_bits")
    )
    hits = _filter_blast_results(
        hits,
        lambda hit: hit.identity_p
        >= core_config.getfloat("blast", "exact_match_min_identity_p"),
    )
    hits = _filter_blast_results(
        hits,
        lambda hit: hit.coverage_p
        >= core_config.getfloat("blast", "exact_match_min_query_cov_p"),
    )
    return [accession for hit in hits for accession in hit.accessions]


def run_blast_similarity(sequence: str,
                         executable_path: str,
                         database_path: str,
                         is_nucleotide: bool,
                         working_directory: str) -> list[tuple[str, float]]:
    """
    Executes similarity search with BLAST
    :param sequence: to be searched
    :param executable_path: path to appropriate (protein or nucleotide) BLAST search engine
    :param database_path: path to BLAST database
    :param is_nucleotide: if True, utilizes nucleotide optimization parameters (instead of protein ones)
    :param working_directory: working directory of the analysis
    :return: accessions of similar proteins with identity percentage
    :raises IdentifierResolutionError: on blast failure (finding no results is NOT a failure)
    """
    hits: list[BlastHit] = _process_blast_results(_run_blast_local(sequence, executable_path,
                                                                   database_path, is_nucleotide, working_directory))

    hits = _filter_blast_results(
        hits,
        lambda hit: hit.evalue
        <= core_config.getfloat("blast", "similarity_match_max_evalue"),
    )

    hits = _filter_blast_results(
        hits,
        lambda hit: hit.bits >= core_config.getint("blast", "similarity_match_min_bits"),
    )

    hits = _filter_blast_results(
        hits,
        lambda hit: hit.identity_p
        >= core_config.getfloat("blast", "similarity_match_min_identity_p"),
    )

    hits = _filter_blast_results(
        hits,
        lambda hit: hit.coverage_p
        >= core_config.getfloat("blast", "similarity_match_min_query_cov_p"),
    )
    return [(accession, hit.identity_p) for hit in hits for accession in hit.accessions]
