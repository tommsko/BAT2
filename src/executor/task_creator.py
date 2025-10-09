import logging
import os.path
from enum import Enum
from rich import print
from .task import Task, TaskList
from .. import core_config
from ..loader.results_file import ResultsFile
from ..utils.logger_utils import create_task_log_output, close_task_log_output
from os import listdir
from os.path import isfile, join
import zipfile
import json

from ..signatures.signature_manager import SignatureGeneratorManager
from ..identifiers.identifier_manager import IdentifierResolutionManager
from ..verifier.summarizer import ResultsSummarizer
from ..verifier.metadata_downloader import MetadataDownloader
from ..loader.simulation import Simulation

logger: logging.Logger = logging.getLogger("base")

class TaskType(Enum):
    ALL = 0
    INPUT = 1
    CREATE = 2
    LOAD = 3
    GENERATE_SIGNATURES = 4
    RESOLVE_IDENTIFIERS = 5
    SUMMARIZE = 6
    METADATA = 7
    UNLOAD = 69


def parse_task_types(tasks: list[str], verbose: bool = True) -> set[TaskType]:
    """
    Parses text-like task types (case-insensitive) into task types
    :param tasks: to parse
    :param verbose: print messages about success/errors
    :return: list of TaskTypes
    """

    res: set[TaskType] = set()
    for task_txt in tasks:
        for task_enum in TaskType:
            if task_enum.name.lower() == task_txt.lower():
                res.add(task_enum)
                break
        else:
            if verbose:
                print(f"[red][b]✘  unable to parse task '{task_txt}', "
                      f"available options: {[x.name for x in TaskType]}[/b][/red]")
            return set()

    if TaskType.ALL in res:
        res = set([x for x in TaskType if x != TaskType.ALL])

    if verbose:
        print(f"[green][b]✔  Following tasks will be executed: "
              f"{[x.name for x in sorted(res, key=lambda x: x.value)]}[/b][/green]")
    return res


def unpack_archives_in_directory(dir_path: str, verbose: bool) -> None:
    """
    Unpacks any .zip files within the directory
    :param dir_path: in which to search for & unpack zip files
    :param verbose: if True, emits messages about unpacking files
    :return: None
    """

    zip_files: list[str] = [os.path.abspath(join(dir_path, f))
                            for f in listdir(dir_path) if isfile(join(dir_path, f)) and f.endswith(".zip")]

    for zip_file in zip_files:
        if verbose:
            print(f"[gray]...  Unpacking '{os.path.basename(zip_file)}' in {os.path.basename(dir_path)}")
        with zipfile.ZipFile(zip_file, 'r') as zip_ref:
            zip_ref.extractall(dir_path)


def save_files_in_directory(path: str) -> None:
    """
    Saves list of all files in directory into the file "all_files.txt"
    :param path: to the directory
    :return: None (creates the file)
    """
    files = [f for f in listdir(path) if isfile(join(path, f))]
    with open(join(path, "all_files.txt"), "w") as f:
        json.dump(files, f)


def create_tasks_for_simulations(simulation_folders: list[str],
                                 task_types: set[TaskType],
                                 verbose: bool = True) -> list[Task]:
    """
    Creates tasks for executor given simulation folders to process and selected tasks
    :param simulation_folders: folders containing simulation data (each to be analyzed)
    :param task_types: to be executed on simulations
    :param verbose: print messages about success/errors
    :return: list of tasks for executor core
    """
    res: list[Task] = []

    for folder in simulation_folders:

        task_log_file_path: str = os.path.join(folder, "tasks.log")
        task_progress_info: dict[str, str] = {"path": task_log_file_path}

        # begin with extracting any zip archives
        unpack_archives_in_directory(folder, verbose)

        save_files_in_directory(folder)

        tasks: list[Task] = []
        def mk_task_list() -> None:
            """
            Appends current tasks to the result as a tasklist
            :return: None
            """
            task_list: TaskList = TaskList(tasks=tasks,
                                           chain_output=True,
                                           timeout_s=core_config.getint("tasks", "global_timeout_s"),
                                           name=f"Process simulation '{folder}'",
                                           task_progress_info=task_progress_info)
            res.append(task_list)

        if TaskType.INPUT in task_types:
            # -> sim_folder
            # <- sim_folder
            tasks.append(Task(executable=lambda x: x,
                              params=[folder],
                              timeout_s=core_config.getint("tasks", "input_timeout_s"),
                              name=f"Input simulation '{folder}'",
                              task_progress_info=task_progress_info))
        else:
            continue  # not even input

        if TaskType.CREATE in task_types:
            # -> results_path, log_handler, sim_path (chained)
            # <- simulation

            def create_sim_task(res, path) -> Simulation:
                log_file_path: str = os.path.join(folder, "results.log")
                log_handler = create_task_log_output(log_file=log_file_path)
                logger.info("Opening task log file '%s'" % log_file_path)
                return Simulation(path, res, log_handler, log_file_path)

            results_file_path: str = os.path.join(folder, "results.json")
            tasks.append(Task(executable=create_sim_task,
                              params=[results_file_path],
                              timeout_s=core_config.getint("tasks", "create_timeout_s"),
                              name=f"Create process for simulation '{folder}'",
                              task_progress_info=task_progress_info))
        else:
            mk_task_list()
            continue

        if TaskType.LOAD in task_types:
            # -> simulation (chained)
            # <- simulation

            def load_sim(sim: Simulation) -> Simulation:
                """
                Loads the simulation and chains it to the output
                :param sim: Simulation to load
                :return: the same simulation as inputted, but loaded
                """
                sim.load_simulation()
                return sim

            tasks.append(
                Task(executable=lambda x: load_sim(x),  # input simulation (chained)
                     params=[],
                     timeout_s=core_config.getint("tasks", "load_timeout_s"),
                     name=f"Load simulation '{folder}'",
                     task_progress_info=task_progress_info))
        else:
            mk_task_list()
            continue

        if TaskType.GENERATE_SIGNATURES in task_types:
            # -> simulation (chained)
            # <- simulation
            tasks.append(
                Task(executable=lambda sim: SignatureGeneratorManager().generate_signatures(sim),
                     params=[],
                     timeout_s=core_config.getint("tasks", "signatures_timeout_s"),
                     name=f"Generate fragment signatures for simulation '{folder}'",
                     task_progress_info=task_progress_info))
        else:
            mk_task_list()
            continue

        if TaskType.RESOLVE_IDENTIFIERS in task_types:
            # -> simulation (chained)
            # <- simulation
            tasks.append(
                Task(executable=lambda sim: IdentifierResolutionManager().resolve_identifiers(sim),
                     params=[],
                     timeout_s=core_config.getint("tasks", "identifiers_timeout_s"),
                     name=f"Resolve signatures to identifiers for simulation '{folder}'",
                     task_progress_info=task_progress_info))
        else:
            mk_task_list()
            continue

        if TaskType.SUMMARIZE in task_types:
            # -> simulation (chained)
            # <- simulation

            tasks.append(
                Task(executable=lambda sim: ResultsSummarizer.summarize(sim),
                     params=[],
                     timeout_s=core_config.getint("tasks", "summarize_timeout_s"),
                     name=f"Summarize results for simulation '{folder}'",
                     task_progress_info=task_progress_info))
        else:
            mk_task_list()
            continue

        if TaskType.METADATA in task_types:
            # -> simulation (chained)
            # <- simulation
            tasks.append(
                Task(executable=lambda sim: MetadataDownloader.fetch_metadata(sim),
                     params=[],
                     timeout_s=core_config.getint("tasks", "metadata_timeout_s"),
                     name=f"Download metadata for results for simulation '{folder}'",
                     task_progress_info=task_progress_info))
        else:
            mk_task_list()
            continue

        ... # any new tasks here

        if TaskType.UNLOAD in task_types:
            # -> simulation (chained)
            # <- simulation

            def unload_sim(sim: Simulation) -> Simulation:
                """
                Unloads the simulation and chains it to the output
                :param sim: Simulation to load
                :return: the same simulation as inputted, but loaded
                """
                sim.unload_simulation()
                if sim.results_file.is_log_handler_attached():
                    logger.info("Closing task log file '%s'" % sim.log_path)
                    close_task_log_output(sim.results_file.get_log_handler())

                errors = []
                with open(sim.log_path, 'r') as fd:
                    for line in fd:
                        if "ERROR" in line or "CRITICAL" in line:
                            errors.append(line)

                error_file_path: str = os.path.join(folder, "errors.log")
                with open(error_file_path, "w") as fd:
                    for line in errors:
                        fd.write(line + "\n")

                return sim

            tasks.append(
                Task(executable=lambda x: unload_sim(x),  # input simulation (chained)
                     params=[],
                     timeout_s=core_config.getint("tasks", "load_timeout_s"),
                     name=f"Unload simulation '{folder}'",
                     task_progress_info=task_progress_info))
        else:
            mk_task_list()
            continue

        mk_task_list()

    if verbose:
        print(f"[green][b]✔  Created {len(res)} tasks to execute")
    return res
