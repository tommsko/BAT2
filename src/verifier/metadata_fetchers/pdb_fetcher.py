import os.path
from typing import Optional, Tuple
import pubchempy as pcp
from chembl_webresource_client.new_client import new_client
import requests
import pyvips
import xml.etree.ElementTree as ET

from .metadata_fetcher import MetadataFetcher
from ...loader.results_file import ResultsFile

class PDBFetcher(MetadataFetcher):
    @staticmethod
    def fetch_metadata(ident: str, results: ResultsFile, fragment_name: str) -> bool:
        """
        Fetches metadata from an identifier.
        :param ident: to be fetched
        :param results: results file to save metadata to
        :param fragment_name: for saving in the correct results section
        :return: True if successful, False otherwise
        """
        ident = ident.upper()
        if '_' in ident:
            ident = ident.split('_')[0]

        results.set_item(f"summary/segments/{fragment_name}/db_crosslink",
                         f"https://www.rcsb.org/structure/{ident}")
        results.set_item(f"summary/segments/{fragment_name}/db_entry_id", ident)
        results.set_item(f"summary/segments/{fragment_name}/db_type", "PDB")

        image_path: str = os.path.join(os.path.dirname(results.path), f"{fragment_name}.jpeg")
        img_data = requests.get(
            f"https://cdn.rcsb.org/images/structures/{ident.lower()}_assembly-1.jpeg").content
        with open(image_path, 'wb') as handler:
            handler.write(img_data)
        results.set_item(f"summary/segments/{fragment_name}/image_path", os.path.basename(image_path))

        full_name = ident  # fallback if API fails
        try:
            url = f"https://www.ebi.ac.uk/proteins/api/proteins/pdb:{ident}"
            r = requests.get(url, headers={"Accept": "application/xml"}, timeout=10)
            r.raise_for_status()

            # Parse XML
            ns = {"u": "https://uniprot.org/uniprot"}
            root = ET.fromstring(r.text)
            fn_elem = root.find(".//u:fullName", ns)
            if fn_elem is not None and fn_elem.text:
                full_name = fn_elem.text.strip()
            else:
                ns = {"u": "http://uniprot.org/uniprot"}
                root = ET.fromstring(r.text)
                fn_elem = root.find(".//u:fullName", ns)
                if fn_elem is not None and fn_elem.text:
                    full_name = fn_elem.text.strip()
        except Exception as e:
            print(f"Warning: could not fetch full name for {ident}: {e}")

        results.set_item(f"summary/segments/{fragment_name}/name", full_name)
        results.set_item(f"summary/segments/{fragment_name}/accession", f"PDB: {ident}")
        return True