from multiprocessing import Process
import logging
import time
from collections import deque
from threading import Thread
from typing import Any, Callable

from ..utils.singleton import Singleton
from .. import core_config

from .task import Task
from .task_creator import TaskType

EMIT_STATS_DELAY_SECONDS = 5

logger: logging.Logger = logging.getLogger("base")


class PseudoProcess(Process):
    """
    A wrapper for "process" to be run in main thread
    """
    def __init__(self, target: Callable[[], Any]):
        """
        Initializes the "pseudo" process
        :param task: executable to be run
        """
        super().__init__()
        self.executable: Callable[[], Any] = target
        self._exitcode: int | None = None

    def start(self) -> None:
        """
        Executes the task, setting internal exitcode
        :return: None
        """
        try:
            self.executable()
            self._exitcode = 0
        except Exception:
            self._exitcode = 1

    def is_alive(self) -> bool:
        """
        The "pseudo" process is never alive, as it runs in the same thread - thus checking means it has finished
        :return: False
        """
        return False

    def join(self, timeout = None) -> None:
        """
        No-op
        :return: None
        """
        ...

    @property
    def exitcode(self) -> int | None:
        return self._exitcode


class ExecutorCore(metaclass=Singleton):
    def __init__(self,
                 terminate_on_finished: bool,
                 block_until_finished: bool,
                 force_single_thread: bool = False) -> None:
        """
        Multiprocess executor manager for tasks
        :param terminate_on_finished: when all tasks are finished, terminate the executor
        :param block_until_finished: run() will block until executor finishes
        :param force_single_thread: force single thread execution (BAT won't spawn any threads)
        """
        self.max_threads: int = core_config.getint("executor", "n_threads")
        self.force_single_thread: bool = force_single_thread

        self._pool_map: dict[Process, Task] = dict()  # mapping active processes to tasks
        self._queue: deque[Task] = deque()  # queue of tasks to be executed

        self.terminate_on_finished: bool = terminate_on_finished
        self.block_until_finished: bool = block_until_finished

        self._active: bool = False
        self._worker_thread: Thread | None = None

        self._worker_sleep_s: float = 1 / core_config.getint("executor", "TPS")
        self._worker_emit_stats: bool = core_config.getboolean("executor", "emit_queue_stats")

    def add_task(self, task: Task) -> str:
        """
        Enqueues task for execution
        :param task: to be executed
        :return: uid of the task
        """
        self._queue.append(task)
        return task.uid

    def get_task_position(self, uid: str) -> int:
        """
        Retrieves task position in the queue
        :param uid: unique identifier of the task
        :return: position in the queue, -1 if not found
        """

        for task in self._pool_map.values():
            if task.uid == uid:
                return 0

        for i in range(len(self._queue)):
            if self._queue[i].uid == uid:
                return i

        return -1

    def run(self) -> None:
        """
        Starts the executor
        :return: None
        """
        self._active = True

        if self.force_single_thread:
            logger.warning("[Executor] Executor will execute all tasks in the main thread (force_single_thread == True)")
            self._worker()
        else:
            self._worker_thread = Thread(target=self._worker)
            self._worker_thread.daemon = True  # make sure thread exits if the app gets killed
            self._worker_thread.start()

            if self.block_until_finished:
                self._worker_thread.join()
                self._worker_thread = None
                self.shutdown()

    def shutdown(self) -> None:
        """
        Terminates the executor
        :return: None
        """
        self._active = False
        self._queue.clear()

        if self._worker_thread is not None:
            self._worker_thread.join()

        for process in self._pool_map.keys():
            if process.is_alive():
                process.kill()
        self._pool_map.clear()

    def is_active(self) -> bool:
        """
        :return: True if executor is still working
        """
        return (len(self._queue) > 0 or len(self._pool_map) > 0) and self._worker_thread.is_alive()

    def _worker_collect(self) -> None:
        """
        Collect already finished tasks
        :return: None
        """
        for process in list(self._pool_map.keys()):
            if process.is_alive():
                continue
            process.join(timeout=1)

            task: Task = self._pool_map[process]
            if process.exitcode == 0:
                logger.debug(f"[Executor] Task {task.uid}{': ' + task.name if task.name else ''}"
                            f" finished successfully!")
                task.dump_task_progress()
            else:
                logger.error(f"[Executor] Task {task.uid}{': ' + task.name if task.name else ''}"
                             f" failed! Check logs")

                for task_type in TaskType:
                    if task_type.name != "ALL":
                        task.task_progress_info[task_type.name] = "FAILED (timeout)"
                task.dump_task_progress()

            del self._pool_map[process]

    def _worker_emit(self) -> None:
        """
        Adds new tasks to process pool while there are tasks in queue and maximum concurrent jobs has not been reached
        :return: None
        """
        while self._queue and len(self._pool_map) < self.max_threads:
            task: Task = self._queue.popleft()
            process: Process = Process(target=task.execute) if not self.force_single_thread \
                                else PseudoProcess(target=task.execute)
            process.start()
            self._pool_map[process] = task
            logger.debug(f"[Executor] Task {task.uid}{': ' + task.name if task.name else ''}"
                         f" is being placed for execution (pid {process.pid})!")

    def _worker(self) -> None:
        """
        Asynchronous worker for enqueueing new tasks
        :return: None
        """

        if self._worker_emit_stats:
            logger.info("[Executor] Executor is starting!")

        last_emitted_msg_time: float = 0

        while self._active:
            if self._worker_emit_stats and time.time() - last_emitted_msg_time > EMIT_STATS_DELAY_SECONDS:
                last_emitted_msg_time = time.time()
                logger.info(f"Executor is running! "
                            f"workers: {self.max_threads}, "
                            f"active: {len(self._pool_map)}, "
                            f"queue: {len(self._queue)}")

            if not self._queue and not self._pool_map and self.terminate_on_finished:
                if self._worker_emit_stats:
                    logger.info("Executor is shutting down! Reason: all jobs finished")
                return

            self._worker_collect()
            self._worker_emit()
            time.sleep(self._worker_sleep_s)
