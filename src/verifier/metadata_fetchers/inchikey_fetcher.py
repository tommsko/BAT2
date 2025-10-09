import os.path
from typing import Optional, Tuple
import pubchempy as pcp
from chembl_webresource_client.new_client import new_client
import requests
import pyvips

from .metadata_fetcher import MetadataFetcher
from ...loader.results_file import ResultsFile

DB_LINK = str
DB_ENTRYID = str
DB_TYPE = str
IMAGE_PATH = str
NAME = str

class InchiKeyFetcher(MetadataFetcher):

    @staticmethod
    def _fetch_pubchem(inchikey: str, img_save_path: str) -> Optional[Tuple[DB_LINK, DB_ENTRYID, DB_TYPE, IMAGE_PATH, NAME]]:
        try:
            compounds = pcp.get_compounds(inchikey, namespace='inchikey')
            if not compounds:
                return None

            compound = compounds[0]
            cid = compound.cid
            name: str = compound.iupac_name or compound.synonyms[0] if compound.synonyms else "N/A"

            img_data = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/image/imagefly.cgi?cid={cid}&width=500&height=500").content
            with open(img_save_path, 'wb') as handler:
                handler.write(img_data)

            return f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}", cid, "PUBCHEM", img_save_path, name
        except Exception:
            return None

    @staticmethod
    def _fetch_chembl(inchikey: str, img_save_path: str) -> Optional[Tuple[DB_LINK, DB_ENTRYID, DB_TYPE, IMAGE_PATH, NAME]]:
        try:
            compound = new_client.molecule.get(inchikey)
            if not compound:
                return None

            chembl_id = compound['molecule_chembl_id']
            name = compound.get('pref_name') or "N/A"

            img_data_raw = requests.get(f"https://www.ebi.ac.uk/chembl/api/data/image/{chembl_id}.svg")
            if img_data_raw.status_code != 200:
                return None
            img = pyvips.Image.svgload_buffer(img_data_raw.content)
            img.write_to_file(img_save_path)
            return f"https://www.ebi.ac.uk/chembl/explore/compound/{chembl_id}", chembl_id, "CHEMBL", img_save_path, name

        except Exception:
            return None

    @staticmethod
    def fetch_metadata(ident: str, results: ResultsFile, fragment_name: str) -> bool:
        """
        Fetches metadata from an identifier.
        :param ident: to be fetched
        :param results: results file to save metadata to
        :param fragment_name: for saving in the correct results section
        :return: True if successful, False otherwise
        """
        image_path: str = os.path.join(os.path.dirname(results.path), f"{fragment_name}.jpg")

        data: Optional[Tuple[DB_LINK, DB_ENTRYID, DB_TYPE, IMAGE_PATH, NAME]] = InchiKeyFetcher._fetch_pubchem(ident, image_path)
        if not data:
            data = InchiKeyFetcher._fetch_chembl(ident, image_path)
        if not data:
            return False

        db_link, db_entryid, db_type, image_path, name = data

        if name == "oxidane":  # bugfix :)
            name = "water"
            db_link = "https://pubchem.ncbi.nlm.nih.gov/compound/962"
            db_entryid = "962"
            db_type = "PUBCHEM"

        results.set_item(f"summary/segments/{fragment_name}/db_crosslink", db_link)
        results.set_item(f"summary/segments/{fragment_name}/db_entry_id", db_entryid)
        results.set_item(f"summary/segments/{fragment_name}/db_type", db_type)
        results.set_item(f"summary/segments/{fragment_name}/image_path", os.path.basename(image_path))
        results.set_item(f"summary/segments/{fragment_name}/name", name)

        results.set_item(f"summary/segments/{fragment_name}/accession", f"{db_type}: {db_entryid}")
        return True