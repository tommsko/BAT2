import logging

import MDAnalysis

from ... import core_config
from ...loader.simulation import Simulation, SegmentSupports
from .signature_generator import SignatureGenerator, SignatureGeneratorType, SignatureGenerationError
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, EditableMol, Mol


class SimpleSmilesGenerator(SignatureGenerator):
    def __init__(self, sanitize: bool = True) -> None:
        """
        Initializes SMILE signature generator
        :param sanitize: of false, generated molecule will NOT be sanitized for chemical sanity
        """
        if sanitize:
            super().__init__("smiles",
                             SignatureGeneratorType.ANALYTICAL,
                             "smiles")
        else:
            super().__init__("smiles_not_sanitized",
                             SignatureGeneratorType.ANALYTICAL,
                             "smiles")

        self._sanitize: bool = sanitize


    def _is_applicable(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> bool:
        """
        Checks whether the fragment can be described by SMILES string
        :param simulation: simulation containing the fragment
        :param fragment: molecule to be described by the SMILES signature
        :return: True if SMILES signature is applicable, False otherwise
        """
        fragment_name: str = fragment.segids[0]

        # we need atomic elements
        if not simulation.get_segment_flag(fragment_name, SegmentSupports.ELEMENTS):
            if not simulation.fix_missing_elements(fragment_name):
                return False

        # we need at least two atoms (aka, a compound) and all of them interconnected
        if len(fragment.atoms.ids) < 2:
            return False

        # proteins being described using SMILES string tremendously slow down the analysis & timeout it
        # we want to prevent this from happening
        if len(fragment.atoms.ids) > core_config.getint("signatures", "smiles_maximum_number_of_atoms"):
            return False

        atom_indices: dict[int, int] = {
            int(atom_idx): i for i, atom_idx in enumerate(fragment.atoms.indices)
        }
        is_atom_bonded: list[bool] = [False for _ in range(len(atom_indices))]
        for atom_idx1, atom_idx2 in fragment.bonds.indices:
            is_atom_bonded[atom_indices[int(atom_idx1)]] = True
            is_atom_bonded[atom_indices[int(atom_idx2)]] = True
        if not all(is_atom_bonded):
            return False

        return True


    def _generate_signature(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> str:
        """
        Generates a SMILES signature for the given fragment
        :param simulation: simulation containing the fragment
        :param fragment: molecule to be described as SMILES string
        :raises SignatureGenerationError: if the generation fails
        :return: SMILES signature
        """

        external_elements: list[str] | None = None
        if any(x == '' for x in fragment.atoms.elements):
            # non-deterministic bug in MDAnalysis when newly set elements 'disappear'
            external_elements = []
            simulation.fix_missing_elements(fragment.segids[0], external_elements)

        atom_idx_to_indices: dict[int, int] = {
            atom_idx: i for i, atom_idx in enumerate(fragment.atoms.indices)
        }

        molecule_build: EditableMol = Chem.EditableMol(Chem.Mol())
        atom_idx_to_rdkit_idx: dict[int, int] = _fragment_to_rdkit_atoms(fragment,
                                                                         atom_idx_to_indices,
                                                                         molecule_build,
                                                                         external_elements)
        _fragment_to_rdkit_bonds(fragment, atom_idx_to_rdkit_idx, molecule_build)

        molecule: Mol = molecule_build.GetMol() # hydrogens are already present from simulation
        _fix_bond_orders_rdkit(molecule)
        molecule = Chem.RemoveHs(molecule, sanitize=self._sanitize)

        if self._sanitize:
            try:
                Chem.SanitizeMol(molecule)
            except Exception as exc:
                raise SignatureGenerationError(f"Failed to sanitize molecule: {exc}") from exc

        return Chem.MolToSmiles(molecule)

class SimpleSmilesGeneratorNoSanitize(SimpleSmilesGenerator):
    def __init__(self) -> None:
        """
        Initializes SMILE signature generator w/o sanitization (wrapper for SimpleSmilesGenerator)
        """
        super().__init__(sanitize=False)

def _fragment_to_rdkit_atoms(
    fragment: MDAnalysis.AtomGroup,
    atom_idx_to_indices: dict[int, int],
    molecule: EditableMol,
    external_elements: list[str] | None,
) -> dict[int, int]:
    """
    Adds all MDAnalysis' fragment atoms to the RDKit's molecule
    :param fragment: AtomGroup representing a molecule
    :param atom_idx_to_indices: mapping of atom indexes to their position in (fragment's) arrays
    :param molecule: RDKit's molecule to add atoms to
    :param external_elements: if provided, overrides list of existing elements from mdanalysis
    :return: mapping of MDAnalysis' atom indexes to RDKit's atom indexes
    """

    idx_atom_to_rdkit: dict[int, int] = {}
    for atom_idx in fragment.atoms.indices:
        if external_elements:
            atom_elem = str(external_elements[atom_idx_to_indices[atom_idx]]).title()
        else:
            atom_elem = str(fragment.atoms.elements[atom_idx_to_indices[atom_idx]]).title()
        atom = Chem.Atom(atom_elem)
        atom.SetNoImplicit(True)
        idx_atom_to_rdkit[atom_idx] = molecule.AddAtom(atom)
    return idx_atom_to_rdkit


def _fragment_to_rdkit_bonds(
    fragment: MDAnalysis.AtomGroup,
    atom_idx_to_rdkit_idx: dict[int, int],
    molecule: EditableMol,
) -> None:
    """
    Adds all MDAnalysis bonds to RDKit's molecule
    :param fragment: AtomGroup representing a molecule
    :param atom_idx_to_rdkit_idx: mapping of MDAnalysis' atom indexes to RDKit's atom indexes
    :param molecule: RDKit's molecule to add bonds to
    :return: None
    """

    for atom_idx1, atom_idx2 in fragment.bonds.indices:
        atom_idx1_rdkit, atom_idx2_rdkit = (
            atom_idx_to_rdkit_idx[atom_idx1],
            atom_idx_to_rdkit_idx[atom_idx2],
        )
        molecule.AddBond(atom_idx1_rdkit, atom_idx2_rdkit, Chem.BondType.SINGLE)

def _fix_bond_orders_rdkit(molecule: Mol,
                           max_attempts: int = 20) -> None:
    """
    Attempts to determine bond orders using RDKit's DetermineBondOrders
    :param max_attempts: maximum attempts to call DetermineBondOrders
    :param molecule: RDKit's molecule
    :return: None (throws on failure)
    """
    charge: int = sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())
    while max_attempts > 0:
        max_attempts -= 1
        try:
            rdDetermineBonds.DetermineBondOrders(
                molecule, allowChargedFragments=True, embedChiral=False, charge=charge
            )
        except ValueError as e:
            print(e)
            charge = int(str(e).split("input (")[1].split(");")[0])
        else:
            return
    raise SignatureGenerationError(
        "unable to determine DetermineBonds() charge parameter within max attempts!"
    )