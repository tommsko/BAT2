import logging
import requests
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from .verification_core import VerificationWorker

logger: logging.Logger = logging.getLogger("verification")

def pubchem_fetch_smiles(cid: str):
    """
    Attempts to fetch SMILES from PubChem.
    :param cid: identifier in PUBCHEM
    :return: SMILES if available, else None
    """
    url: str = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/TXT"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        return r.text.strip()
    except Exception as e:
        logger.warning(f"PubChem {cid} fetch SMILES failed: {e}")
        return None

def chembl_fetch_smiles(chembl_id: str):
    """
    Attempts to fetch smiles given CHEMBL ID.
    :param chembl_id: identifier in CHEMBL (without initial CHEMBL string)
    :return: SMILES if available, else None
    """

    url: str = f"https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL{chembl_id}.json"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        data = r.json()
        return data.get("molecule_structures", {}).get("canonical_smiles")
    except Exception as e:
        logger.warning(f"CHEMBL {chembl_id} fetch SMILES failed: {e}")
        return None

def chebi_fetch_smiles(chebi_id: str):
    """
    Attempts to fetch SMILES given CHEBI database.
    :param chebi_id: identifier in CHEBI (without initial CHEBI string)
    :return: SMILES if available, else None
    """

    chebi_num: str = chebi_id.replace("CHEBI:", "")
    url: str = f"https://www.ebi.ac.uk/chebi/ws/rest/chemicalStructure/{chebi_num}"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        data = r.text
        if "<smiles>" in data:
            smiles = data.split("<smiles>")[1].split("</smiles>")[0]
            return smiles.strip()
        return None
    except Exception as e:
        logger.warning(f"ChEBI {chebi_id} fetch SMILES failed: {e}")
        return None

def tanimoto_similarity(smiles1: str, smiles2: str):
    """
    Compares two SMILES using tanimoto similarity.
    :param smiles1: molecule 1
    :param smiles2: molecule 2
    :return: tanimoto similarity
    """
    if not smiles1 or not smiles2:
        return 0.0
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if mol1 is None or mol2 is None:
        return 0.0
    fp1 = FingerprintMols.FingerprintMol(mol1)
    fp2 = FingerprintMols.FingerprintMol(mol2)
    return DataStructs.FingerprintSimilarity(fp1, fp2) * 100.0

class CompoundVerificationWorker(VerificationWorker):
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

        if ident_type == "PUBCHEM_CID":
            seq = pubchem_fetch_smiles(ident)
        elif ident_type == "CHEMBL_ID":
            seq = chembl_fetch_smiles(ident)
        elif ident_type == "CHEBI_ID":
            seq = chebi_fetch_smiles(ident)
        else:
            logger.error(f"Unknown compound ident_type: {ident_type}")
            seq = None

        if seq is None or not seq:
            return None

        return tanimoto_similarity(signature, seq)
