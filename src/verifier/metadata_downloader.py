import logging
from . import VerificationError
from ..loader.results_file import ResultsFile
import os
from os import listdir
from os.path import isfile, join
from .metadata_fetchers.export_simulation_snapshot import generate_snapshot
from .metadata_fetchers.atom_fetcher import AtomFetcher
from .metadata_fetchers.inchikey_fetcher import InchiKeyFetcher
from .metadata_fetchers.pdb_fetcher import PDBFetcher
from ..loader.simulation import Simulation

logger: logging.Logger = logging.getLogger("verification")

metadata_resolvers: dict[str, type | None] = {
    "BLAST_PDB_PROT": PDBFetcher,
    "BLAST_PDB_NT": PDBFetcher,
    "BLAST_UNIPROT": None,
    "PDB_STRUCT": PDBFetcher,
    "ALPHAFIND_STRUCT": PDBFetcher,
    "CHEMBL_STRUCT": InchiKeyFetcher,
    "PUBCHEM_STRUCT": InchiKeyFetcher,
    "ATOM": AtomFetcher,
    "NAME_CCD": InchiKeyFetcher,
    "NAME_MAPS": InchiKeyFetcher,
    "NAME_MANUAL": InchiKeyFetcher,
}

class MetadataDownloader:
    # downloads basic metadata for found hits

    @staticmethod
    def _mk_placeholders(results: ResultsFile, fragment_name: str) -> None:
        """
        Creates placeholders (None values) for to-be-fetched fields
        :param results: results file associated with the simulation
        :param fragment_name: to create placeholders for
        :return: None
        """
        results.set_item(f"summary/segments/{fragment_name}/db_crosslink", None)
        results.set_item(f"summary/segments/{fragment_name}/db_entry_id", None)
        results.set_item(f"summary/segments/{fragment_name}/db_type", None)

        results.set_item(f"summary/segments/{fragment_name}/image_path", None)

        results.set_item(f"summary/segments/{fragment_name}/name", None)
        results.set_item(f"summary/segments/{fragment_name}/accession", None)

    @staticmethod
    def _export_simulation_snapshot(results: ResultsFile, sim: Simulation) -> bool:
        sim_files: list[str] = [f for f in listdir(sim.simulation_directory)
                                if isfile(join(sim.simulation_directory, f))]
        tpr_path = [path for path in sim_files if path.endswith(".tpr")]
        if not tpr_path:
            return False
        tpr_path = os.path.join(sim.simulation_directory, tpr_path[0])

        gro_path = [path for path in sim_files if path.endswith(".gro")]
        if not gro_path:
            gro_path = None
        else:
            gro_path = os.path.join(sim.simulation_directory, gro_path[0])

        image_path: str = os.path.join(os.path.dirname(results.path), "SIMULATION.png")
        generate_snapshot(tpr_path, gro_path, image_path)

        results.set_item(f"summary/segments/SIMULATION/db_crosslink",None)
        results.set_item(f"summary/segments/SIMULATION/db_entry_id", None)
        results.set_item(f"summary/segments/SIMULATION/db_type", None)
        results.set_item(f"summary/segments/SIMULATION/image_path", "SIMULATION.png")
        results.set_item(f"summary/segments/SIMULATION/name", "Automatic annotation of the simulation")
        results.set_item(f"summary/segments/SIMULATION/identifier", "BAT2")
        results.set_item(f"summary/segments/SIMULATION/ident", os.path.basename(tpr_path))
        results.set_item(f"summary/segments/SIMULATION/confidence", None)

        return True

    @staticmethod
    def fetch_metadata(sim: Simulation) -> Simulation:
        """
        Fetches basic metadata for summarized segments
        :param results_path: results file associated with the simulation
        :param sim: unused, part of task chaining
        :return: None
        """

        # results: ResultsFile = sim.get_results_handler()
        #
        # segments: list[str] = results.get_item("summary/segment_list")
        # for segment in segments:
        #
        #     if segment == "SIMULATION":
        #         continue
        #
        #     MetadataDownloader._mk_placeholders(results, segment)
        #
        #     identifier: str | None = results.get_item(f"summary/segments/{segment}/identifier")
        #     if identifier is None:
        #         continue
        #     ident: str = results.get_item(f"summary/segments/{segment}/ident")
        #
        #     if metadata_resolvers.get(identifier) is None:
        #         logger.error(f"Undefined metadata fetcher for identifier '{identifier}'")
        #         continue
        #
        #     metadata_resolvers[identifier].fetch_metadata(ident, results, segment)

        # replaced with mol*
        # MetadataDownloader._export_simulation_snapshot(results, sim)

        return sim

