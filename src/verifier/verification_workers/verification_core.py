from enum import Enum

from src.loader.results_file import ResultsFile
from src.loader.simulation import Simulation
from src import core_config
from tqdm import tqdm

class VerificationKind(Enum):
    PROTEIN = 1
    NUCLEIC = 2
    COMPOUND = 3
    WATER = 4
    ATOM = 5
    NONE = 6


class VerificationWorker:
    def verify_identification(self, ident: str, ident_type: str, signature: str) -> float | None:
        """
        Verification of identification using reverse database lookup
        :param ident: identifier to be verified
        :param ident_type: type of identifier
        :param signature: signature to verify identification against
        :return: new confidence given verification, or None if verification failed
        """
        raise NotImplementedError  # to be implemented by derived classes


class Verifier:
    def __init__(self, simulation: Simulation, results: ResultsFile,
                 fragment_name: str, verification_kind: VerificationKind, worker: VerificationWorker | None) -> None:
        """
        Initializes results verifier
        :param simulation: being verified
        :param results: results file object
        :param fragment_name: fragment name
        :param verification_kind: type of verification to be executed
        :param worker: verification worker to process identifiers and signatures into observed confidence
        """
        self.simulation: Simulation = simulation
        self.results: ResultsFile = results
        self.fragment_name: str = fragment_name
        self.verification_kind: VerificationKind = verification_kind
        self.worker: VerificationWorker | None = worker

    def get_signature(self) -> str | None:
        """
        Returns a signature of the fragment
        :return: FASTA residue sequence for protein, SMILES string for compounds, None on failure
        """
        all_signatures: dict[str, str] | None = self.results.get_item(f"signatures/{self.fragment_name}/generated",
                                                                default=None)
        if all_signatures is None:
            return None

        for k, v in all_signatures.items():
            if self.verification_kind == VerificationKind.PROTEIN and k in ('FASTA_MAP_PROTEIN', 'FASTA_RESNAME_PROTEIN'):
                return v
            if self.verification_kind == VerificationKind.NUCLEIC and k in ('FASTA_MAP_NUCLEIC', 'FASTA_RESNAME_NUCLEIC'):
                return v
            if self.verification_kind == VerificationKind.COMPOUND and k in ('SMILES', 'SMILES_NOT_SANITIZED'):
                return v
        return None

    def get_identifiers(self) -> list[tuple[str, str, float]]:
        """
        Returns available identifiers for the fragment
        :return: list of (identifier, identifier_type, confidence), sorted by highest confidence first
        """

        results: list[tuple[str, str, float]] = []
        for signature_gen, data in self.results.get_item(f"identifiers/{self.fragment_name}").items():
            for identifier, data_ident in data.items():
                if identifier == "attempted":
                    continue

                ident_type = data_ident["hit_type"]
                idents = data_ident["hits"]
                for ident, confidence in idents:
                    results.append((ident, ident_type, float(confidence) if isinstance(confidence, str) else confidence))
        return sorted(results, key=lambda x: x[2], reverse=True)

    def calculate_identifier_confidence(self) -> list[tuple[str, str, float, float | None]]:
        """
        Calculates verified confidence for all identifiers, if possible
        :return: list (entry for each identifier) of (identifier, identifier_type, confidence, observed_confidence),
                    sorted by observed confidence
        """

        identifiers: list[tuple[str, str, float]] = self.get_identifiers()

        if self.verification_kind == VerificationKind.WATER:
            return [("962", "PUBCHEM_CID", 100.0, 100.0)]

        if self.verification_kind == VerificationKind.ATOM:
            return [(self.get_identifiers()[0][0], "ELEMENT", 100.0, 100.0)]

        if self.verification_kind == VerificationKind.NONE:
            identifiers.sort(key=lambda x: x[2], reverse=True)
            return [(i, it, sc, None) for i, it, sc in identifiers]

        results: list[tuple[str, str, float, float | None]] = []

        to_verify: int = core_config.getint("verification", "verify_top_n_results")
        successfully_verified: int = 0
        for identifier, identifier_type, confidence in tqdm(identifiers, desc=f"Verifying {self.fragment_name} identifiers..."):
            observed_confidence: float | None = self.worker.verify_identification(identifier, identifier_type, self.get_signature())
            results.append((identifier, identifier_type, confidence, observed_confidence if observed_confidence is not None else -1))

            if observed_confidence is not None:
                successfully_verified += 1

            if successfully_verified > to_verify:
                break

        results.sort(key=lambda x: x[3], reverse=True)
        return [(i, it, sc, oc if oc != -1 else None) for i, it, sc, oc in results]

