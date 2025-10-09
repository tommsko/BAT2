import logging
from collections import deque
from typing import Callable

import MDAnalysis

from ... import core_config
from ...loader.simulation import Simulation, SegmentSupports
from .signature_generator import SignatureGenerator, SignatureGeneratorType, SignatureGenerationError
from ...constants import KNOWN_RESIDUE_CODES, KNOWN_NUCLEIC_ACIDS_CODES


class MapFastaGenerator(SignatureGenerator):
    def __init__(self, is_protein: bool, is_nucleotide: bool) -> None:
        """
        Initializes FASTA signature generator
        :param is_protein: initialize as protein descriptor
        :param is_nucleotide: initialize as nucleotide descriptor
        """

        assert is_protein or is_nucleotide, "One mode of operation is mandatory"
        assert not (is_protein and is_nucleotide), "Only one mode of operation is allowed"
        self._is_protein: bool = is_protein
        self._is_nucleotide: bool = is_nucleotide

        super().__init__("fasta_map_protein" if is_protein else "fasta_map_nucleic",
                         SignatureGeneratorType.ANALYTICAL,
                         "residue_sequence_protein" if is_protein else "residue_sequence_nucleotide")


    def _is_applicable(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> bool:
        """
        Checks whether the fragment can be described by an amino/nucleic acid residue sequence
        :param simulation: simulation containing the fragment
        :param fragment: molecule to be described by a residue sequence
        :return: True if fasta signature is applicable, False otherwise
        """
        fragment_name: str = fragment.segids[0]

        # we need correctly named residues (fixing not yet supported)
        if self._is_protein:
            return (simulation.get_segment_flag(fragment_name, SegmentSupports.AMINO_ACID_RESIDUE_NAMES)
                    and simulation.get_segment_flag(fragment_name, SegmentSupports.MULTI_ATOM_RESIDUES))
        else: # nucleic
            return (simulation.get_segment_flag(fragment_name, SegmentSupports.NUCLEIC_ACID_RESIDUE_NAMES)
                    and not simulation.get_segment_flag(fragment_name, SegmentSupports.MULTI_ATOM_RESIDUES))


    def _generate_signature(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> str:
        """
        Generates an amino/nucleic acid residue sequence signature for the given fragment (using map aproach)
        :param simulation: simulation containing the fragment
        :param fragment: molecule to be described by a residue sequence
        :raises SignatureGenerationError: if the generation fails
        :return: residue sequence (fasta-like)
        """

        return fragment_to_protseq_map_based(fragment)

class MapFastaGeneratorProtein(MapFastaGenerator):
    def __init__(self) -> None:
        super().__init__(is_protein=True, is_nucleotide=False)

class MapFastaGeneratorNucleotide(MapFastaGenerator):
    def __init__(self) -> None:
        super().__init__(is_protein=False, is_nucleotide=True)

class ResnameFastaGenerator(SignatureGenerator):
    def __init__(self, is_protein: bool, is_nucleotide: bool) -> None:
        """
        Initializes FASTA signature generator
        :param is_protein: initialize as protein descriptor
        :param is_nucleotide: initialize as nucleotide descriptor
        """

        assert is_protein or is_nucleotide, "One mode of operation is mandatory"
        assert not (is_protein and is_nucleotide), "Only one mode of operation is allowed"
        self._is_protein: bool = is_protein
        self._is_nucleotide: bool = is_nucleotide

        super().__init__("fasta_resname_protein" if is_protein else "fasta_resname_nucleic",
                         SignatureGeneratorType.ANALYTICAL,
                         "residue_sequence_protein" if is_protein else "residue_sequence_nucleotide")

    def _is_applicable(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> bool:
        """
        Checks whether the fragment can be described by an amino/nucleic acid residue sequence
        :param simulation: simulation containing the fragment
        :param fragment: molecule to be described by a residue sequence
        :return: True if fasta signature is applicable, False otherwise
        """
        fragment_name: str = fragment.segids[0]

        # we need correctly named residues (fixing not yet supported)
        if self._is_protein:
            return (simulation.get_segment_flag(fragment_name, SegmentSupports.AMINO_ACID_RESIDUE_NAMES)
                    and simulation.get_segment_flag(fragment_name, SegmentSupports.MULTI_ATOM_RESIDUES))
        else: # nucleic
            return (simulation.get_segment_flag(fragment_name, SegmentSupports.NUCLEIC_ACID_RESIDUE_NAMES)
                    and not simulation.get_segment_flag(fragment_name, SegmentSupports.MULTI_ATOM_RESIDUES))


    def _generate_signature(self, simulation: Simulation, fragment: MDAnalysis.AtomGroup) -> str:
        """
        Generates an amino/nucleic acid residue sequence signature for the given fragment (using map aproach)
        :param simulation: simulation containing the fragment
        :param fragment: molecule to be described by a residue sequence
        :raises SignatureGenerationError: if the generation fails
        :return: residue sequence (fasta-like)
        """

        return fragment_to_protseq_resname_based(fragment)

class ResnameFastaGeneratorProtein(ResnameFastaGenerator):
    def __init__(self) -> None:
        super().__init__(is_protein=True, is_nucleotide=False)

class ResnameFastaGeneratorNucleotide(ResnameFastaGenerator):
    def __init__(self) -> None:
        super().__init__(is_protein=False, is_nucleotide=True)

def resname_to_one_letter(resname: str) -> str:
    """
    Transforms residue name (amino or nucleic acid) to one-letter representation
    :param resname: name of the residue representing either amino or nucleic acid
    :return: one-letter residue code
    """
    if len(resname) == 1:
        return resname.upper()  # presumes it is already a code, i.e. nucleic acids

    if resname.upper() in KNOWN_NUCLEIC_ACIDS_CODES:
        return KNOWN_NUCLEIC_ACIDS_CODES[resname.upper()]

    return KNOWN_RESIDUE_CODES[resname.upper()]

def fragment_to_protseq_resname_based(fragment: MDAnalysis.AtomGroup) -> str:
    """
    Assumes the order of residue indices is the same as the order of protein residues and extracts sequence based
    on this assumption
    :param fragment: molecule to be described by an amino/nucleic acid sequence
    :raises SignatureGenerationError: if the fragment description by residue sequence fails
    :return: residue sequence (fasta-like)
    """

    try:
        residue_atom_indices: list[int] = [  # residues consist of multiple atoms (with same resname), we want only one
            indices[0] for indices in fragment.residues.indices
        ]

        atom_indices_to_index: dict[int, int] = {  # mapping of atom indices (mdanalysis) to their positions in arrays
            indice: i for i, indice in enumerate(fragment.atoms.indices)
        }

        residue_atom_resnames: list[str] = [  # transform residue names into a one-letter codes for fasta-like output
            resname_to_one_letter(fragment.resnames[atom_indices_to_index[i]].upper())
            for i in residue_atom_indices
        ]

        residues_and_indices: list[tuple[int, str]] = list(
            zip(residue_atom_indices, residue_atom_resnames)
        )
        residues_and_indices.sort(key=lambda x: x[0])  # sort by indices
        return "".join(x[1] for x in residues_and_indices)

    except Exception as exc:
        raise SignatureGenerationError("protein fingerprint resolution (resname-based)") from exc

def fragment_to_protseq_map_based(fragment: MDAnalysis.AtomGroup) -> str:
    """
    Tries to translate fragment into amino/nucleic acid residue sequence utilizing graph approach (safer)
    :param fragment: molecule to be described by an amino/nucleic acid sequence
    :raises SignatureGenerationError: if the fragment description by residue sequence fails
    :return: residue sequence (fasta-like)
    """
    prot_sign_allow_partial_seq: bool = (
        core_config.getboolean('signatures', 'fasta_map_determine_sequence_partial_mode'))
    prot_sign_strict_mode: bool = (
        core_config.getboolean('signatures', 'fasta_map_determine_sequence_strict_mode'))

    (
        atom_indices_to_index,
        atom_adjacency_matrix,
        residue_indices_to_code,
        residue_indices_to_index,
        residue_adjacency_matrix,
    ) = _fragment_to_protein_mk_adjacency_matrices(fragment)

    _fragment_to_protein_sequence_verify_sequentiality(
        residue_adjacency_matrix, atom_adjacency_matrix
    )

    methionine_indices, non_methionine_indices = _split_residue_indices_methionine(
        residue_indices_to_code
    )

    methionine_sequences: list[str] = _find_longest_sequences_bfs(
        methionine_indices, residue_adjacency_matrix, residue_indices_to_index, residue_indices_to_code
    )
    if methionine_sequences:
        if sel := _filter_sequences(
            methionine_sequences, require=lambda seq: len(seq) == len(residue_indices_to_code)
        ):
            sel.sort(key=lambda seq: len(seq), reverse=True)
            return sel[0]

    logging.getLogger("signature").debug(
        "... unable to determine protein sequence given exact & strict mode. "
        "If allowed, continuing exact & non-strict mode"
    )
    if prot_sign_strict_mode:
        cant_determine_err(prot_sign_strict_mode, prot_sign_allow_partial_seq)

    non_methionine_sequences: list[str] = _find_longest_sequences_bfs(
        non_methionine_indices, residue_adjacency_matrix, residue_indices_to_index, residue_indices_to_code
    )
    if sel := _filter_sequences(
        non_methionine_sequences, require=lambda seq: len(seq) == len(residue_indices_to_code)
    ):
        sel.sort(key=lambda seq: len(seq), reverse=True)
        return sel[0]

    logging.getLogger("signature").debug(
        "... unable to determine protein sequence given exact & non-strict mode. "
        "If allowed, continuing partial & non-strict mode"
    )
    if not prot_sign_allow_partial_seq:
        cant_determine_err(prot_sign_strict_mode, prot_sign_allow_partial_seq)

    all_sequences = methionine_sequences + non_methionine_sequences

    all_sequences.sort(key=lambda seq: len(seq), reverse=True)
    if all_sequences:
        return all_sequences[0]

    cant_determine_err(prot_sign_strict_mode, prot_sign_allow_partial_seq)
    assert False  # unreachable

################################################################
# helper functions for map-based approach

AdjacencyMatrix = list[list[bool]]


def _fragment_extract_residues(
    fragment: MDAnalysis.AtomGroup,
) -> tuple[dict[int, str], dict[int, int]]:
    """
    Extracts amino/nucleic acid residues from a fragment
    :param fragment: AtomGroup representing a molecule to be described by amino/nucleic acid sequence
    :return: mapping from residue indices (MDAnalysis) to 1-letter residues' code
            and mapping of residue index to residue indices
    """

    residue_indices_to_code: dict[int, str] = {}
    residue_indices_to_index: dict[int, int] = {}
    res_i: int = 0
    for res_name, res_indice in zip(fragment.resnames, fragment.resindices):  # residue indices repeat themselves
        res_name = resname_to_one_letter(res_name)                            # one residue has multiple atoms
        if res_indice not in residue_indices_to_code:                         # (which share resname and resindice)
            residue_indices_to_code[res_indice] = res_name
            residue_indices_to_index[res_indice] = res_i
            res_i += 1
    return residue_indices_to_code, residue_indices_to_index


def _mk_adjacency_matrix(size: int) -> AdjacencyMatrix:
    """
    Creates adjacency matrix
    :param size: size of the (square) matrix
    :return: adjacency matrix
    """
    return [[False for _ in range(size)] for _ in range(size)]


def _populate_adjacency_matrix_inplace(
    fragment: MDAnalysis.AtomGroup,
    atom_indices_to_index: dict[int, int],
    residue_indices_to_index: dict[int, int],
    adjacency_matrix: AdjacencyMatrix,
    is_residue_adjacency: bool,
) -> None:
    """
    Populates adjacency matrix for given bonds in the fragment
    :param fragment: AtomGroup representing a molecule to be described by amino/nucleic acid sequence
    :param atom_indices_to_index: mapping from atom indices to their positions in arrays
    :param residue_indices_to_index: mapping from residue indices to their positions in arrays
    :param adjacency_matrix: result of _mk_adjacency_matrix(), will be populated
    :param is_residue_adjacency: if True, populates it with residue adjacency, otherwise atom adjacency
    :return: None
    """

    for atom_idx1, atom_idx2 in fragment.bonds.indices:
        atom_i1, atom_i2 = atom_indices_to_index[atom_idx1], atom_indices_to_index[atom_idx2]
        atom_resid_i1, atom_resid_i2 = (  # less like residue index and more like order
            residue_indices_to_index[fragment.atoms.resindices[atom_i1]],
            residue_indices_to_index[fragment.atoms.resindices[atom_i2]],
        )

        node_i1, node_i2 = (
            (atom_resid_i1, atom_resid_i2)
            if is_residue_adjacency
            else (atom_i1, atom_i2)
        )
        adjacency_matrix[node_i1][node_i2] = adjacency_matrix[node_i2][node_i1] = True


def _adjacency_matrix_get_neighbors(
    node_from: int, adjacency_matrix: AdjacencyMatrix
) -> set[int]:
    """
    Finds all linked nodes (neighbors) to a nodes
    :param node_from: node's index
    :param adjacency_matrix: adjacency matrix between the nodes
    :return: set of all neighbors, except queried node
    """
    result: set[int] = set()
    for neighbor in range(len(adjacency_matrix)):
        if node_from == neighbor:
            continue
        if (
            adjacency_matrix[node_from][neighbor]
            or adjacency_matrix[neighbor][node_from]
        ):
            result.add(neighbor)
    return result


def _bfs(
    node_start: int,
    adjacency_matrix: AdjacencyMatrix,
    reuse_visited: list[bool] | None = None,
) -> tuple[list[bool], list[int], list[int], int]:
    """
    Generic BFS implementation
    :param node_start: node index from which to start BFS
    :param adjacency_matrix: adjacency matrix of nodes
    :param reuse_visited: if provided, re-uses visited matrix
    :return: list of visited nodes (possibly reused),
             list of distances to start_node,
             list of predecessors given BFS search,
             number of nodes visited during BFS search
    """

    visited: list[bool] = (
        [False for _ in range(len(adjacency_matrix))]
        if reuse_visited is None
        else reuse_visited
    )
    distances: list[int] = [-1 for _ in range(len(adjacency_matrix))]
    prev: list[int] = [-1 for _ in range(len(adjacency_matrix))]
    nodes_visited: int = 0

    queue: deque[int] = deque()
    queue.append(node_start)
    visited[node_start] = True
    distances[node_start] = 0
    prev[node_start] = node_start

    while queue:
        node_from: int = queue.popleft()
        nodes_visited += 1
        for node_to in _adjacency_matrix_get_neighbors(node_from, adjacency_matrix):
            if visited[node_to]:
                continue
            visited[node_to] = True
            distances[node_to] = distances[node_from] + 1
            prev[node_to] = node_from
            queue.append(node_to)

    return visited, distances, prev, nodes_visited


def _discontinuous_segment_cnt(adjacency_matrix: AdjacencyMatrix) -> int:
    """
    Finds how many weakly-linked components there in adjacency matrix
    :param adjacency_matrix: to be analysed
    :return: number of weakly-linked components
    """
    visited = [False for _ in range(len(adjacency_matrix))]
    components: int = 0
    for i_from in range(len(adjacency_matrix)):
        if visited[i_from]:
            continue
        _bfs(i_from, adjacency_matrix, visited)
        components += 1
    return components


def _degrees_count(adjacency_matrix: AdjacencyMatrix) -> list[int]:
    """
    Calculates sum of in and out degree count for each node in adjacency matrix
    :param adjacency_matrix: to be analysed
    :return: degrees for each node
    """
    degrees: list[int] = [0 for _ in range(len(adjacency_matrix))]
    for node_from in range(len(adjacency_matrix)):
        for node_to in range(len(adjacency_matrix)):
            if node_from != node_to and adjacency_matrix[node_from][node_to]:
                degrees[node_from] += 1  # symmetric
    return degrees


def _bfs_extract_longest_paths(
    node_start: int,
    distances: list[int],
    prev: list[int],
    visited: list[bool],
    max_length: int = 10000,
) -> list[list[int]]:
    """
    Extracts the path(s) between node_start and most distant node(s) utilizing results of BFS calculation before
    :param node_start: the residue *INDICE, NOT IDX* from which BFS was run
    :param distances: distances list from BFS run
    :param prev: previous list from BFS run
    :param visited: visited list from BFS run
    :param max_length: maximum length of the path (to throw on cycles)
    :return: path from node_start to node with the biggest distance, or more if there are multiple nodes with such distance
    """
    biggest_distance: int = max(distances)
    nodes_to: list[int] = [
        node_i
        for node_i, node_distance in enumerate(distances)
        if node_distance == biggest_distance
    ]
    results: list[list[int]] = []

    for node_to in nodes_to:
        path: list[int] = []
        while node_to != node_start:
            assert visited[node_to], "corrupted BFS calculation"
            path.append(node_to)
            node_to = prev[node_to]
            if len(path) > max_length:
                raise SignatureGenerationError(
                    "finding residue sequence: "
                    f"residues' path is too long or cyclic (max_length = {max_length})",
                )
        path.append(node_start)
        results.append(list(reversed(path)))
    return results


def _split_residue_indices_methionine(
    residues: dict[int, str]
) -> tuple[set[int], set[int]]:
    """
    Splits residues' indices to two groups, methionine and non-methionine ones
    :param residues: result from _fragment_extract_residues(...)
    :return: set of methionine indices, set of non-methionine indices
    """
    is_methionine: set[int] = set()
    not_methionine: set[int] = set()

    for residue_id, residue_code in residues.items():
        if residue_code == "M":
            is_methionine.add(residue_id)
        else:
            not_methionine.add(residue_id)
    return is_methionine, not_methionine


def _translate_residue_indices_to_sequence(
    indices_path: list[int], residue_indices_to_index: dict[int, int], residue_indices_to_code: dict[int, str]
) -> str:
    """
    Translates residue indices path from BFS to sequence of 1-code amino/nucleic acids
    :param indices_path: path of residue indices
    :param residue_indices_to_index: mapping of residue indices to their index
    :param residue_indices_to_code: mapping of residues indices to their one-letter code
    :return: sequence string
    """
    res_i_to_idx: dict[int, int] = {
        resid_i: resid_idx for resid_idx, resid_i in residue_indices_to_index.items()
    }
    return "".join(residue_indices_to_code[res_i_to_idx[x]] for x in indices_path)


def _find_longest_sequences_bfs(
    from_indices: set[int],
    residue_adjacency_matrix: AdjacencyMatrix,
    residue_indices_to_index: dict[int, int],
    residue_indices_to_code: dict[int, str],
) -> list[str]:
    """
    Finds longest sequences given residues and their adjacency matrix and filters them
    :param from_indices: set of residue IDs from which longest sequence is being searched
    :param residue_adjacency_matrix: adjacency matrix of residues
    :param residue_indices_to_index: mapping from residue indices to index
    :param residue_indices_to_code: mapping of residue indices to their one-letter codes
    :return: list of sequences found
    """
    sequences: list[str] = []
    for residue_idx in from_indices:
        node_from: int = residue_indices_to_index[residue_idx]
        vis, dist, prev, _ = _bfs(node_from, residue_adjacency_matrix)
        paths = _bfs_extract_longest_paths(node_from, dist, prev, vis)
        for path in paths:
            sequences.append(
                _translate_residue_indices_to_sequence(path, residue_indices_to_index, residue_indices_to_code)
            )
    return sequences


def _filter_sequences(
    sequences: list[str], require: Callable[[str], bool]
) -> list[str]:
    """
    Simple wrapper to filter a list of strings
    :param sequences: to be filtered
    :param require: lambda that must hold true to include sequence in the result
    :return: list of sequences where require(...) returns True
    """
    return [seq for seq in sequences if require(seq)]


def _fragment_to_protein_sequence_verify_sequentiality(
    residue_adjacency_matrix: AdjacencyMatrix, atom_adjacency_matrix: AdjacencyMatrix
) -> None:
    """
    Verifies that residues and atoms form one component, and that residues form one sequential chain
    :param residue_adjacency_matrix: adjacency matrix of residues
    :param atom_adjacency_matrix: adjacency matrix of atoms
    :return: None
    :raises SignatureGenerationError: of there are multiple components or degrees mismatch
    """

    if (residue_components := _discontinuous_segment_cnt(residue_adjacency_matrix)) != 1:
        raise SignatureGenerationError(
            "extracting residue chain: "
            f"found {residue_components} disjoint residue components, "
            f"which is more than 1 supported!",
        )
    if (atoms_components := _discontinuous_segment_cnt(atom_adjacency_matrix)) != 1:
        raise SignatureGenerationError(
            "extracting residue chain: "
            f"found {atoms_components} disjoint atom components, "
            f"which is more than 1 supported! (although there's only one residue component)",
        )

    residue_degrees: list[int] = _degrees_count(residue_adjacency_matrix)
    if any(d >= 3 or d == 0 for d in residue_degrees):
        raise SignatureGenerationError(
            "extracting residue chain: found a residue with degree > 2, sequence is non-linear!"
        )
    if residue_degrees.count(1) >= 3:
        raise SignatureGenerationError(
            "extracting residue chain: found more than two residues with degree == 1, sequence is non-linear!",
        )


def _mk_residue_atom_adjacency_matrix(
    fragment: MDAnalysis.AtomGroup,
    residue_indices_to_code: dict[int, str],
    atom_indices_to_index: dict[int, int],
    residue_indices_to_index: dict[int, int],
) -> tuple[AdjacencyMatrix, AdjacencyMatrix]:
    """
    Creates residue and atom adjacency matrix given
    :param fragment: AtomGroup representing a molecule to be described by amino/nucleic acid sequence
    :param residue_indices_to_code: mapping of residue indices to their one-letter codes
    :param atom_indices_to_index: mapping of atom indices to their positions in arrays
    :param residue_indices_to_index: mapping of residue ids to residue indices
    :return: residue and atom adjacency matrix
    """

    residue_adjacency_matrix: AdjacencyMatrix = _mk_adjacency_matrix(len(residue_indices_to_code))
    _populate_adjacency_matrix_inplace(
        fragment,
        atom_indices_to_index,
        residue_indices_to_index,
        residue_adjacency_matrix,
        is_residue_adjacency=True,
    )

    atom_adjacency_matrix: AdjacencyMatrix = _mk_adjacency_matrix(len(atom_indices_to_index))
    _populate_adjacency_matrix_inplace(
        fragment,
        atom_indices_to_index,
        residue_indices_to_index,
        atom_adjacency_matrix,
        is_residue_adjacency=False,
    )
    return residue_adjacency_matrix, atom_adjacency_matrix


def cant_determine_err(strict: bool, partial: bool) -> None:
    """
    Just raises an exception
    :param strict: residue sequence determination is in strict mode
    :param partial: protein sequence determination is in partial mode
    """
    raise SignatureGenerationError(
        "map-based residue sequence signature generation: "
        f"Unable to determine the sequence! (modes: exact=True, strict={strict}, "
        f"partial={partial})",
    )


def _fragment_to_protein_mk_adjacency_matrices(fragment: MDAnalysis.AtomGroup) -> ...:
    """
    Creates indices and adjacency matrices for residues and atoms given
    :param fragment: AtomGroup representing a molecule (protein specifically)
    :return: atom_indices - mapping of atom idx to atom indices
             atom_adjacency_matrix - adjacency matrix of atoms, indices-based
             residues - mapping of residue idx to residue code
             residue_indices - mapping of residue idx to residue indices
             residue_adjacency_matrix - adjacency matrix of residues, indices-based
    """

    atom_indices_to_index: dict[int, int] = {
        indice: index for index, indice in enumerate(fragment.atoms.indices)
    }

    residue_indices_to_code: dict[int, str]
    residue_indices_to_index: dict[int, int]
    residue_indices_to_code, residue_indices_to_index = _fragment_extract_residues(fragment)

    residue_adjacency_matrix: AdjacencyMatrix
    atom_adjacency_matrix: AdjacencyMatrix

    residue_adjacency_matrix, atom_adjacency_matrix = _mk_residue_atom_adjacency_matrix(
        fragment, residue_indices_to_code, atom_indices_to_index, residue_indices_to_index
    )
    return (
        atom_indices_to_index,
        atom_adjacency_matrix,
        residue_indices_to_code,
        residue_indices_to_index,
        residue_adjacency_matrix,
    )
