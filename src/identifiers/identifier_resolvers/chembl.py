
from chembl_webresource_client.settings import Settings

Settings.Instance().TIMEOUT = 10
Settings.Instance().CACHING = False
Settings.Instance().TOTAL_RETRIES = 10
from chembl_webresource_client.new_client import new_client
from .identifier_resolver import IdentifierResolver, IdentifierResolutionError
from ... import core_config


class CHEMBLResolver(IdentifierResolver):
    """
    CHEMBLResolver tries to resolve smiles identifiers using 2D structural CHEMBL database search
    """

    def __init__(self) -> None:
        """
        Initializes CHEMBLResolver for given signatures
        """
        super().__init__("chembl_struct", "CHEMBL_ID",["smiles"])


    def _is_applicable(self, signature: str) -> bool:
        """
        Checks whether signature is SMILES and applicable for selected CHEMBL database
        :param signature: textual signature to be identified
        :return: True if identifier CHEMBL is applicable for signature, False otherwise
        """
        return True  # TODO: presupposing that output of SMILES signature generator is good.
                     # No further checks are needed


    def _resolve_identity(self, signature: str, working_directory: str) -> list[str]:
        """
        Implementation for resolving signature on identity level
        :param signature: textual signature to be resolved
        :param working_directory: unused
        :raises IdentifierResolutionError: if resolution fails
        :return: list of identifiers exactly matched to the signature
        """
        return chembl_fetch_similarity_match(
                signature, threshold=100
            )


    def _resolve_similarity(self, signature: str, working_directory: str) -> list[tuple[str, float]]:
        """
        Implementation for resolving signature on similarity level
        :param signature: textual signature to be resolved
        :param working_directory: unused
        :raises IdentifierResolutionError: if resolution fails (no results are okay)
        :return: list of identifiers with associated similarity scores (0-1)
        """
        return chembl_fetch_similarity_match(signature,
                                             threshold=core_config.getint("chembl",
                                                                          "similarity_threshold"),
                                             include_similarity=True
            )

def chembl_fetch_similarity_match(smiles: str, threshold: int, include_similarity: bool = False) -> list[str | tuple[str, float]]:
    """
    Searches CHEMBL's databases for similarity matches of SMILES
    :param smiles: identifier to be searched
    :param threshold: of similarity to report results (100 == identity)
    :param include_similarity: if True, similarity percentages will be included in the result
    :raises IdentifierResolutionError: on CHEMBL search failure
    :return: list of all InChI Keys of molecules similar to smiles query given threshold
    """

    # https://notebooks.gesis.org/binder/jupyter/user/chembl-chembl_webresource_client-5neb7ifo/notebooks/demo_wrc.ipynb

    try:
        results: list[str | tuple[str, float]] = []
        similarity = new_client.similarity
        records = similarity.filter(smiles=smiles, similarity=threshold).only(
            ["molecule_chembl_id", "similarity"]
        )

        for record in records:
            if (
                "similarity" not in record
                or "molecule_chembl_id" not in record
                or "CHEMBL" not in record["molecule_chembl_id"]
            ):
                continue
            similarity: float = record["similarity"]
            chembl_id: str = record["molecule_chembl_id"].split("CHEMBL")[1]
            if include_similarity:
                results.append((chembl_id, similarity))
            else:
                results.append(chembl_id)

        return results
    except Exception as exc:
        raise IdentifierResolutionError(
            f"[Identifier] CHEMBL searching failed: {exc}"
        )
