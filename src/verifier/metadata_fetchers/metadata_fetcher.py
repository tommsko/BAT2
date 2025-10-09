from ...loader.results_file import ResultsFile


class MetadataFetcher:
    """
    Abstract class for metadata fetchers
    """

    @staticmethod
    def fetch_metadata(ident: str, results: ResultsFile, fragment_name: str) -> bool:
        """
        Fetches metadata from an identifier. To be implemented by specific fetchers
        :param ident: to be fetched
        :param results: results file to save metadata to
        :param fragment_name: for saving in the correct results section
        :return: True if successful, False otherwise
        """
        raise NotImplementedError