import logging
import os
import warnings
from collections import defaultdict
from io import StringIO

import MDAnalysis
from Bio.PDB import Structure, Model, Chain, Residue, Atom, PDBIO, PDBParser
import numpy as np
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from MDAnalysis import AtomGroup

from ...loader.simulation import Simulation, SegmentSupports
from .signature_generator import SignatureGenerator, SignatureGeneratorType, SignatureGenerationError
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, EditableMol, Mol
from ...utils.file_utils import mk_random_filename
from ... import core_config


class PDBGenerator(SignatureGenerator):
    def __init__(self, backbone_only: bool = False, sanitize: bool = True, is_protein: bool = False, is_nucleotide: bool = False) -> None:
        """
        Initializes PDB signature generator
        :param backbone_only: only backbone (alfa-carbons) will be exported
        :param sanitize: of false, generated pdb will NOT be sanitized
        :param is_protein: initialize as protein descriptor
        :param is_nucleotide: initialize as nucleotide descriptor
        """

        assert is_protein or is_nucleotide, "One mode of operation is mandatory"
        assert not (is_protein and is_nucleotide), "Only one mode of operation is allowed"
        self._is_protein: bool = is_protein
        self._is_nucleotide: bool = is_nucleotide

        if backbone_only:
            if sanitize:
                super().__init__("pdb_backbone_protein" if is_protein else "pdb_backbone_nucleic",
                                 SignatureGeneratorType.ANALYTICAL,
                                 "pdb_protein" if is_protein else "pdb_nucleic")
            else:
                super().__init__("pdb_backbone_not_sanitized_protein" if is_protein else "pdb_backbone_not_sanitized_nucleic",
                                 SignatureGeneratorType.ANALYTICAL,
                                 "pdb_protein" if is_protein else "pdb_nucleic")
        else:
            if sanitize:
                super().__init__("pdb_protein" if is_protein else "pdb_nucleic",
                                 SignatureGeneratorType.ANALYTICAL,
                                 "pdb_protein" if is_protein else "pdb_nucleic")
            else:
                super().__init__("pdb_not_sanitized_protein" if is_protein else "pdb_not_sanitized_nucleic",
                                 SignatureGeneratorType.ANALYTICAL,
                                 "pdb_protein" if is_protein else "pdb_nucleic")

        self._backbone_only: bool = backbone_only
        self._sanitize: bool = sanitize


    def _is_applicable(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> bool:
        """
        Checks whether the fragment can be described by PDB file
        :param simulation: simulation containing the fragment
        :param fragment: molecule to be described by the PDF file
        :return: True if PDB signature is applicable, False otherwise
        """
        fragment_name: str = fragment.segids[0]

        # we need atomic positions
        if not simulation.get_segment_flag(fragment_name, SegmentSupports.POSITIONS):
            return False

        # we need to be working with amino/nucleic acid residues
        if self._is_protein:
            if (not simulation.get_segment_flag(fragment_name, SegmentSupports.AMINO_ACID_RESIDUE_NAMES)
                    or not simulation.get_segment_flag(fragment_name, SegmentSupports.MULTI_ATOM_RESIDUES)):
                return False
        else:
            if (not simulation.get_segment_flag(fragment_name, SegmentSupports.NUCLEIC_ACID_RESIDUE_NAMES)
                    or simulation.get_segment_flag(fragment_name, SegmentSupports.MULTI_ATOM_RESIDUES)):
                return False

        # if we are not generating alpha-carbon backbone, we need all elements
        if not self._backbone_only and not simulation.get_segment_flag(fragment_name, SegmentSupports.ELEMENTS):
            if not simulation.fix_missing_elements(fragment_name):
                return False

        # we need unique atomic names in each residue
        if not simulation.get_segment_flag(fragment_name, SegmentSupports.NAMES_UNIQUE):
            if not simulation.fix_missing_elements(fragment_name):
                return False

        return True



    def _generate_signature(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> str:
        """
        Generates a PDB signature for the given fragment
        :param simulation: simulation containing the fragment
        :param fragment: molecule to be described by PDB file
        :raises SignatureGenerationError: if the generation fails
        :return: PDB signature
        """

        pdb: str = fragment_to_pdb(fragment, self._backbone_only)
        if (self._sanitize
                and not verify_pdb(pdb, core_config.get("miscellaneous", "temporary_directory"))):
            raise SignatureGenerationError("Generated PDB couldn't be parsed")
        return pdb

class PDBGeneratorSanitizeProtein(PDBGenerator):
    def __init__(self) -> None:
        """
        Initializes PDB signature generator w/ sanitization (wrapper for PDBGenerator)
        """
        super().__init__(backbone_only=False, sanitize=True, is_protein=True)

class PDBGeneratorNoSanitizeProtein(PDBGenerator):
    def __init__(self) -> None:
        """
        Initializes PDB signature generator w/o sanitization (wrapper for PDBGenerator)
        """
        super().__init__(backbone_only=False, sanitize=False, is_protein=True)

class PDBBackboneGeneratorSanitizeProtein(PDBGenerator):
    def __init__(self) -> None:
        """
        Initializes PDB backbone signature generator w/ sanitization (wrapper for PDBGenerator)
        """
        super().__init__(backbone_only=True, sanitize=True, is_protein=True)

class PDBBackboneGeneratorNoSanitizeProtein(PDBGenerator):
    def __init__(self) -> None:
        """
        Initializes PDB backbone signature generator w/o sanitization (wrapper for PDBGenerator)
        """
        super().__init__(backbone_only=True, sanitize=False, is_protein=True)

class PDBGeneratorSanitizeNucleic(PDBGenerator):
    def __init__(self) -> None:
        """
        Initializes PDB signature generator w/ sanitization (wrapper for PDBGenerator)
        """
        super().__init__(backbone_only=False, sanitize=True, is_nucleotide=True)

class PDBGeneratorNoSanitizeNucleic(PDBGenerator):
    def __init__(self) -> None:
        """
        Initializes PDB signature generator w/o sanitization (wrapper for PDBGenerator)
        """
        super().__init__(backbone_only=False, sanitize=False, is_nucleotide=True)

class PDBBackboneGeneratorSanitizeNucleic(PDBGenerator):
    def __init__(self) -> None:
        """
        Initializes PDB backbone signature generator w/ sanitization (wrapper for PDBGenerator)
        """
        super().__init__(backbone_only=True, sanitize=True, is_nucleotide=True)

class PDBBackboneGeneratorNoSanitizeNucleic(PDBGenerator):
    def __init__(self) -> None:
        """
        Initializes PDB backbone signature generator w/o sanitization (wrapper for PDBGenerator)
        """
        super().__init__(backbone_only=True, sanitize=False, is_nucleotide=True)

def fragment_to_pdb(fragment: MDAnalysis.AtomGroup, only_backbone: bool) -> str:
    """
    Attempts to convert fragment into PDB file using various techniques
    :param fragment: molecule to be described by PDB file
    :param only_backbone: export only alpha-carbons (one atom per residue)
    :raises SignatureGenerationError: if the generation fails
    :return: PDB string
    """

    try:
        return fragment_to_pdb_biopython(fragment, only_backbone)
    except Exception:
        ...

    try:
        return fragment_to_pdb_direct(fragment, only_backbone, include_bonds=True)
    except Exception:
        ...

    try:
        return fragment_to_pdb_direct(fragment, only_backbone, include_bonds=False)
    except Exception:
        raise SignatureGenerationError("Unable to generate PDB signature")

def verify_pdb(pdb: str, temp_directory: str) -> bool:
    """
    Verifies whether PDB file can be parsed
    :param pdb: pdb files content
    :param temp_directory: where to temporarily save pdb file
    :return: True if pdb content is parse-able
    """

    path: str = mk_random_filename("pdb", temp_directory)
    with open(path, 'w') as fd:
        fd.write(pdb)

    success: bool = True
    try:
        warnings.filterwarnings("error", category=PDBConstructionWarning)
        p = PDBParser()
        p.get_structure("export", path)
        warnings.resetwarnings()
    except Exception:
        success = False

    os.unlink(path)
    return success


def fragment_to_pdb_biopython(fragment: MDAnalysis.AtomGroup, only_backbone: bool) -> str:
    """
    Generates (contents of) PDB file from a fragment using biopython
    :param fragment: molecule to be described by PDB file
    :param only_backbone: export only alpha-carbons (one atom per residue)
    :raises SignatureGenerationError: if the generation fails
    :return: PDB string
    """

    atom_ids_to_indices: dict[int, int] = {indice: i for i, indice in enumerate(fragment.atoms.indices)}

    structure: Structure = Structure.Structure("export")
    model: Model = Model.Model(0)
    chain: Chain = Chain.Chain("A")

    for resid, resname in zip(fragment.residues.resnums, fragment.residues.resnames):
        atoms: AtomGroup = fragment.select_atoms(f"resid {resid}")
        residue: Residue = Residue.Residue((" ", resid, " "), resname, "")

        if only_backbone:
            mean_position = np.array(atoms.atoms.positions, dtype=float)
            atom = Atom.Atom(atoms.atoms.names[0],
                             np.array(mean_position.mean(axis=0)),
                             1.0,
                             1.0,
                             " ",
                             "CA",
                             atom_ids_to_indices[atoms.atoms.ids[0]],
                             "C")
            residue.add(atom)
        else:
            for i in range(len(atoms)):
                atom = Atom.Atom(atoms.names[i],
                                 np.array(atoms.positions[i]),
                                 1.0,
                                 1.0,
                                 " ",
                                 atoms.names[i],
                                 atom_ids_to_indices[atoms.ids[0]],
                                 str(atoms.elements[i]).upper())
                residue.add(atom)

        chain.add(residue)

    model.add(chain)
    structure.add(model)

    io: PDBIO = PDBIO()
    io.set_structure(structure)

    pdb_string_io = StringIO()
    io.save(pdb_string_io)
    pdb: str = pdb_string_io.getvalue()
    pdb_string_io.close()
    return pdb

def _pdb_atom_line(
    atom_idx: int,
    atom_name: str,
    atom_elem: str,
    residue_name: str,
    residue_idx: int,
    pos_x: np.float32,
    pos_y: np.float32,
    pos_z: np.float32,
) -> str:
    """
    Creates well-formatted ATOM line of .pdb file
    :param atom_idx: atom ID
    :param atom_name: atom name (must be unique!)
    :param atom_elem: atom element
    :param residue_name: residue name
    :param residue_idx: residue ID
    :param pos_x: atom X position
    :param pos_y: atom Y position
    :param pos_z: atom Z position
    :return: ATOM line with values provided
    """
    # https://cupnet.net/pdb-format/
    return "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(
        "ATOM",
        atom_idx,
        atom_name,
        "",
        residue_name.upper(),
        "A",
        residue_idx,
        "",
        pos_x,
        pos_y,
        pos_z,
        1,
        0,
        atom_elem,
        "",
    )


def _pdb_bonds_line(atom_idx: int, other_atoms: set[int]) -> str:
    """
    Creates well-formatted CONECT line of .pdb file
    :param atom_idx: atom ID of which bonds are described
    :param other_atoms: IDs of other atoms bonded to this atom
    :return: CONECT line with values provided
    """
    # https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html
    conect_fmts = [
        "{:6s}{:5d}{:5d}",
        "{:6s}{:5d}{:5d}{:5d}",
        "{:6s}{:5d}{:5d}{:5d}{:5d}",
        "{:6s}{:5d}{:5d}{:5d}{:5d}{:5d}",
    ]
    return conect_fmts[len(other_atoms) - 1].format("CONECT", atom_idx, *other_atoms)


def fragment_to_pdb_direct(
    fragment: MDAnalysis.AtomGroup,
    only_backbone: bool,
    include_bonds: bool,
) -> str:
    """
    Converts a fragment into PDB file
    :param fragment: AtomGroup representing a protein
    :param only_backbone: export only alpha-carbons (one atom per residue)
    :param include_bonds: if True, bonds (CONECT lines) will be included
    :raises SignatureGenerationError: if the generation fails
    :return: PDB string
    """
    processed_residues: set[int] = set()

    fd: StringIO = StringIO()

    fd.write("HEADER    export\n")

    for atom_i in range(fragment.n_atoms):
        if only_backbone:
            if fragment.atoms.resindices[atom_i] in processed_residues:
                continue
            else:
                processed_residues.add(fragment.atoms.resindices[atom_i])
        fd.write(
            _pdb_atom_line(
                fragment.atoms.ids[atom_i],
                fragment.atoms.names[atom_i],
                fragment.atoms.elements[atom_i] if not only_backbone else "C",
                fragment.atoms.resnames[atom_i],
                fragment.atoms.resindices[atom_i],
                fragment.atoms.positions[atom_i][0],
                fragment.atoms.positions[atom_i][1],
                fragment.atoms.positions[atom_i][2],
            )
            + "\n"
        )

    if include_bonds:
        conns: dict[int, set] = defaultdict(set)
        for atom_idx1, atom_idx2 in fragment.atoms.bonds.indices:
            conns[atom_idx1].add(atom_idx2)
            conns[atom_idx2].add(atom_idx1)
        for atom_idx in conns:
            fd.write(_pdb_bonds_line(atom_idx, conns[atom_idx]) + "\n")

        fd.write("END\n")

    pdb: str = fd.getvalue()
    fd.close()
    return pdb




