import logging
import os.path
import signal
import dill # type: ignore
import pickle
import json
from typing import Callable, Any, Union

from ..utils.string_utils import random_uid
from ..utils.time_utils import current_sec_time

logger: logging.Logger = logging.getLogger("base")

class Task:
    def __init__(self,
                 executable: Union[Callable[[Any], Any], Callable[[], Any]],
                 params: list[Any],
                 timeout_s: int | None = None,
                 name: str = "",
                 task_progress_info: dict[str, str] | None = None) -> None:
        """
        Creates a generic task to be executed
        :param executable: callable to the task executable, it is expected to throw exception on failure
        :param params: parameters to be passed to the executable
        :param timeout_s: if not None, timeouts **THE PROCESS** after <<timeout_s>> seconds
        :param name: name of the task
        """
        self.uid: str = random_uid()
        self.name: str = name

        self._executable: bytes = dill.dumps(executable)
        self._params: list[Any] = params
        self._timeout_s: int | None = timeout_s

        self.finished_successfully: bool = False

        assert task_progress_info is not None, "task_progress_info is required"
        self.task_progress_info: dict[str, str] = task_progress_info

    def execute(self) -> Any:
        """
        Runs the executable in task, updates self.finished_successfully appropriately
        :return: The output of the task or None on failure
        """

        logger.debug(f"[Task {self.uid}{': ' + self.name if self.name else ''}] starting...")
        started_at_sec = current_sec_time()

        if self._timeout_s is not None:
            signal.alarm(self._timeout_s)

        try:
            ret: Any = dill.loads(self._executable)(*self._params)
            self.finished_successfully = True
            if self._timeout_s is not None: signal.alarm(0)

            runtime_sec = round(current_sec_time() - started_at_sec, 2)
            self.task_progress_info[self.name] = f"OK ({runtime_sec} s / {self._timeout_s})"
            return ret
        except Exception as exc:
            logger.error(f"[Task {self.uid}{': ' + self.name if self.name else ''}] failed! {exc}")
            logger.exception(exc)  # TODO: remove
            if self._timeout_s is not None: signal.alarm(0)

            runtime_sec = round(current_sec_time() - started_at_sec, 2)
            self.task_progress_info[self.name] = f"FAILED ({runtime_sec} s / {self._timeout_s}): {exc}"
            return None

    def append_to_params(self, what: Any) -> None:
        """
        Appends a parameter to the parameters passed to the executable
        :param what: parameter to append
        :return: None
        """
        self._params.append(what)

    def dump_task_progress(self) -> None:
        if "path" not in self.task_progress_info:
            return

        path = self.task_progress_info["path"]
        #del self.task_progress_info["path"]

        tasks_to_save = self.task_progress_info.copy()

        if os.path.exists(path):
            with open(path, "r") as fd:
                data =  json.load(fd)
                for k, v in data.items():
                    if k not in tasks_to_save:
                        tasks_to_save[k] = v

        with open(path, "w") as fd:
            json.dump(tasks_to_save, fd)


class TaskList(Task):
    def __init__(self,
                 tasks: list[Task],
                 chain_output: bool = True,
                 timeout_s: int | None = None,
                 name: str = "",
                 task_progress_info: dict[str, str] | None = None) -> None:
        """
        Creates a chained list of tasks to be executed in order
        :param tasks: list of tasks to be executed, in order
        :param chain_output: if True, output of previous task is appended to the next one (as the last parameter)
        :param timeout_s: if not None, timeouts **THE PROCESS** after <<timeout_s>> seconds
        :param name: name of the TaskList
        """
        super().__init__(lambda: True,
                         [],
                         timeout_s,
                         name,
                         task_progress_info)

        self._tasks: list[Task] = tasks
        self._chain_output: bool = chain_output
        self._timeout_sec: int | None = timeout_s

    def execute(self) -> Any:
        """
        Executes the tasks from TaskList in order. If any task fails, successors are not attempted.
        Updates self.finished_successfully appropriately
        :return: The output of the task or None on failure
        """

        logger.debug(f"[TaskList {self.uid}{': ' + self.name if self.name else ''}] starting...")

        if self._timeout_s is not None:
            signal.alarm(self._timeout_s)

        last_output: Any = None

        for task in self._tasks:
            if self._chain_output and last_output is not None:
                task.append_to_params(last_output)

            last_output = task.execute()
            self.dump_task_progress()

            if not task.finished_successfully:
                logger.error(
                    f"[TaskList {self.uid}{' ' + self.name if self.name else ''}] failed! Subtask failed!")
                if self._timeout_s is not None: signal.alarm(0)
                return None

        self.finished_successfully = True
        if self._timeout_s is not None: signal.alarm(0)
        return last_output
