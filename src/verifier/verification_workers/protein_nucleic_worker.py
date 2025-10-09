import logging

import requests
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from .verification_core import VerificationWorker

logger: logging.Logger = logging.getLogger("verification")

def pdb_fetch_sequence(pdb_id: str) -> str | None:
    """
    Attempts to fetch protein sequence given a PDB ID.
    :param pdb_id: referring to protein/nucleic acid
    :return: FASTA sequence, or None on failure
    """
    if '_' in pdb_id: pdb_id = pdb_id.split('_', maxsplit=1)[0]

    url: str = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id.lower()}/1"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        data = r.json()
        return data['entity_poly']['pdbx_seq_one_letter_code_can'].replace("\n", "")
    except Exception as e:
        logger.warning(f"PDB {pdb_id} fetch FASTA failed: {e}")
        return None

def uniprot_fetch_sequence(uniprot_id: str):
    """
    Fetches sequence given a UniProt ID.
    :param uniprot_id: referring to protein/nucleic acid
    :return: FASTA sequence, or None on failure
    """
    url: str = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        lines = r.text.splitlines()
        return ''.join(line.strip() for line in lines if not line.startswith('>'))
    except Exception as e:
        logger.warning(f"UniProt {uniprot_id} fetch FASTA failed: {e}")
        return None

def alphafold_fetch_sequence(alphafold_id: str):
    """
    Attempts to fetch protein sequence given a AlphaFold ID.
    :param alphafold_id: referring to protein
    :return: FASTA sequence, or None on failure
    """
    try:
        uniprot_id: str = alphafold_id.split('-')[1]
        return uniprot_fetch_sequence(uniprot_id)
    except Exception as e:
        logger.warning(f"AlphaFold {alphafold_id} fetch FASTA failed: {e}")
        return None

def seq_identity(seq1: str, seq2: str):
    """
    Compares two sequences using global alignment
    :param seq1: A
    :param seq2: B
    :return: identity percentage
    """
    if not seq1 or not seq2:
        return 0.0
    alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
    aln = alignments[0]
    matches = sum(a == b for a, b in zip(aln.seqA, aln.seqB))
    identity = matches / max(len(seq1), len(seq2))
    return identity * 100.0

class ProteinNucleicVerificationWorker(VerificationWorker):
    def verify_identification(self, ident: str, ident_type: str, signature: str) -> float | None:
        """
        Verification of identification using reverse database lookup
        :param ident: identifier to be verified
        :param ident_type: type of identifier
        :param signature: signature to verify identification against
        :return: new confidence given verification, or None if verification failed
        """
        if signature is None or not signature:
            return None

        if ident_type == "PDB_ID":
            seq = pdb_fetch_sequence(ident)
        elif ident_type == "UNIPROT_ID":
            seq = uniprot_fetch_sequence(ident)
        elif ident_type == "ALPHAFOLD_ID":
            seq = alphafold_fetch_sequence(ident)
        else:
            logger.warning(f"Unknown protein ident_type: {ident_type}")
            seq = None

        if seq is None or not seq:
            return None

        return seq_identity(signature, seq)
