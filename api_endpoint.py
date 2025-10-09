import os
from rich import print

from src import core_config
from src.utils.logger_utils import init_loggers
from src.executor.task_creator import TaskType, parse_task_types, create_tasks_for_simulations
from src.executor.task import Task
from src.executor.core import ExecutorCore


def run_analyze(directory: str, task_uuid: str) -> str | None:
    """
    Runs full BAT analysis (--run) mode
    :param directory: of the files
    :param task_uuid: task uuid
    :return: path to the resulting .json file
    """
    init_loggers(verbose=True, very_verbose=True)
    print("[cyan][b]-> Entering run mode[/b][/cyan]")
    print(f"[b]About to analyse: '{directory}'[/b]")

    tasks_to_exec: set[TaskType] = parse_task_types(["ALL"], verbose=True)
    tasks: list[Task] = create_tasks_for_simulations([directory], tasks_to_exec)
    print(f"[green][b]✔ {len(tasks)} task(s) created[/b][/green]")

    if not core_config.has_section("TMP"):
        core_config.add_section("TMP")
    core_config.set("TMP", "CURR_UUID", task_uuid)

    executor: ExecutorCore = ExecutorCore(terminate_on_finished=True,
                                          block_until_finished=True,
                                          force_single_thread=False)
    for task in tasks:
        executor.add_task(task)

    print("[cyan][b]-> Launching synchronous executor (will block until finished)[/b][/cyan]")
    executor.run()
    print(f"[green][b]✔ executor has exitted[/b][/green]")


    print(f"[b]Collecting results from: '{directory}'[/b]")
    results_path: str = os.path.join(directory, "results.json")
    if not os.path.exists(results_path):
        print("[red][b]✘ the job failed[/b][/red]")
        return None

    return os.path.abspath(results_path)


