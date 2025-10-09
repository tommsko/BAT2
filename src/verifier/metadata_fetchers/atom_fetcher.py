import os.path

from .metadata_fetcher import MetadataFetcher
from ...loader.results_file import ResultsFile
from ...constants import ATOMIC_ELEMENT_TO_PERIODICTABLE_INDICE, ATOMIC_ELEMENT_TO_NAME

from PIL import Image, ImageDraw, ImageFont

class AtomFetcher(MetadataFetcher):

    @staticmethod
    def create_atomic_element_image(element: str, path: str) -> None:
        """
        Creates simple image for atomic element
        :param element: short name of the element
        :param path: to save image to
        :return: None
        """

        width, height = 300, 300
        image = Image.new("RGB", (width, height), "white")

        draw = ImageDraw.Draw(image)
        font =ImageFont.load_default(80)

        bbox = draw.textbbox((0, 0), element, font=font)
        text_width = bbox[2] - bbox[0]
        text_height = bbox[3] - bbox[1]

        x = (width - text_width) / 2
        y = (height - text_height) / 2

        draw.text((x, y), element, fill="black", font=font)
        image.save(path)

    @staticmethod
    def fetch_metadata(ident: str, results: ResultsFile, fragment_name: str) -> bool:
        """
        Fetches metadata from an identifier.
        :param ident: to be fetched
        :param results: results file to save metadata to
        :param fragment_name: for saving in the correct results section
        :return: True if successful, False otherwise
        """
        element_indice: int = ATOMIC_ELEMENT_TO_PERIODICTABLE_INDICE[ident]

        results.set_item(f"summary/segments/{fragment_name}/db_crosslink",
                         f"https://pubchem.ncbi.nlm.nih.gov/element/{element_indice}")
        results.set_item(f"summary/segments/{fragment_name}/db_entry_id", element_indice)
        results.set_item(f"summary/segments/{fragment_name}/db_type", "ATOM")

        image_path: str = os.path.join(os.path.dirname(results.path), f"{fragment_name}.jpg")
        AtomFetcher.create_atomic_element_image(ident, image_path)
        results.set_item(f"summary/segments/{fragment_name}/image_path", os.path.basename(image_path))
        results.set_item(f"summary/segments/{fragment_name}/name", ATOMIC_ELEMENT_TO_NAME[ident])
        results.set_item(f"summary/segments/{fragment_name}/accession", ATOMIC_ELEMENT_TO_NAME[ident])

        return True