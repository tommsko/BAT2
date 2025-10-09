import logging

from .identifier_resolver import IdentifierResolver, IdentifierResolutionError
from ... import core_config

from SPARQLWrapper import SPARQLWrapper, JSON

log = logging.getLogger("identifier")

class IDSMResolver(IdentifierResolver):

    """
    IDSMResolver tries to resolve smiles identifiers using 2D structural CHEMBL/PubChem database search
    """

    def __init__(self, use_pubchem: bool = False, use_chembl: bool = False, use_chebi: bool = False) -> None:
        """
        Initializes IDSMResolver for given signatures
        """

        assert use_chembl or use_pubchem or use_chebi, "One mode of operation is necessary"
        assert not (use_chembl and use_pubchem), "Only one mode of operation is allowed"

        self._use_pubchem: bool = use_pubchem
        self._use_chembl: bool = use_chembl
        self._use_chebi: bool = use_chebi

        super().__init__("idsm_pubchem" if use_pubchem else "idsm_chembl" if use_chembl else "idsm_chebi",
                         "PUBCHEM_CID" if use_pubchem else "CHEMBL_ID" if use_chembl else "CHEBI_ID",
                         ["smiles"])


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
        return idsm_fetch_similarity_match(
                signature, threshold=100, include_similarity=False,
                as_pubchem=self._use_pubchem, as_chembl=self._use_chembl, as_chebi=self._use_chebi
            )


    def _resolve_similarity(self, signature: str, working_directory: str) -> list[tuple[str, float]]:
        """
        Implementation for resolving signature on similarity level
        :param signature: textual signature to be resolved
        :param working_directory: unused
        :raises IdentifierResolutionError: if resolution fails (no results are okay)
        :return: list of identifiers with associated similarity scores (0-1)
        """

        return idsm_fetch_similarity_match(
                signature,
                threshold=core_config.getint("idsm","similarity_threshold"),
                include_similarity=True,
                as_pubchem=self._use_pubchem, as_chembl=self._use_chembl, as_chebi=self._use_chebi
            )

class IDSMResolverPubchem(IDSMResolver):
    def __init__(self):
        super().__init__(use_pubchem=True)

class IDSMResolverChembl(IDSMResolver):
    def __init__(self):
        super().__init__(use_chembl=True)

class IDSMResolverChebi(IDSMResolver):
    def __init__(self):
        super().__init__(use_chebi=True)

def idsm_fetch_similarity_match(smiles: str, threshold: int | None, include_similarity: bool = False,
                                as_pubchem: bool = False, as_chembl = False, as_chebi = False) -> list[str | tuple[str, float]]:
    """
    Searches IDSM's databases for similarity matches of SMILES in chembl/pubchem/chebi databases
    :param smiles: identifier to be searched
    :param threshold: of similarity to report results (100 == identity)
    :param include_similarity: if True, similarity percentages will be included in the result
    :raises IdentifierResolutionError: on IDSM search failure
    :return: list of all database identifiers (PubChem CID's or CHEMBL ID's) of molecules similar to smiles query given threshold
    """

    cutoff: str = "0.8" if threshold is None else str(round(threshold/100, 1))

    assert as_pubchem or as_chembl or as_chebi, "One mode of operation is mandatory"

    # https://notebooks.gesis.org/binder/jupyter/user/chembl-chembl_webresource_client-5neb7ifo/notebooks/demo_wrc.ipynb

    try:
        endpoint = "https://idsm.elixir-czech.cz/sparql/endpoint/idsm"
        query = '''
        PREFIX sachem: <http://bioinfo.uochb.cas.cz/rdf/v1.0/sachem#>
        PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

        SELECT * WHERE {
          [ sachem:compound ?COMPOUND;
            sachem:score ?SCORE ] sachem:similaritySearch [
                  sachem:query "''' + smiles + '''";
                  sachem:cutoff "''' + cutoff + '''"^^xsd:double ].
        }
            '''

        sparql = SPARQLWrapper(endpoint)
        sparql.setQuery(query)
        sparql.setReturnFormat(JSON)

        results = sparql.query().convert()

        ret = []
        for result in results["results"]["bindings"]:
            id = similarity = None
            if as_pubchem and 'pubchem' in result["COMPOUND"]["value"]:
                id = int(result["COMPOUND"]["value"].split("CID")[1])
                similarity = float(result["SCORE"]["value"]) * 100
            elif as_chembl and 'chembl' in result["COMPOUND"]["value"]:
                id = int(result["COMPOUND"]["value"].split("CHEMBL")[1])
                similarity = float(result["SCORE"]["value"]) * 100
            elif as_chebi and 'CHEBI' in result["COMPOUND"]["value"]:
                id = int(result["COMPOUND"]["value"].split("CHEBI_")[1])
                similarity = float(result["SCORE"]["value"]) * 100

            if id is not None and threshold is not None and similarity < threshold:
                continue

            if id is not None:
                ret.append((id, similarity) if include_similarity else id)

        return ret
    except Exception as exc:
        logging.error("IDSM search failed!")
        logging.exception(exc)

        raise IdentifierResolutionError(
            f"[Identifier] IDSM searching failed: {exc}"
        )
