
import time

from ...src.executor.task import Task, TaskList
from ...src.executor.core import ExecutorCore
from .. import delete_file, check_file, write_file, delete_directory, make_directory, TEST_DIR


# test enqueueing + position
def test_executor_enqueue():
    exe = ExecutorCore(terminate_on_finished=True, block_until_finished=False)

    task1: Task = Task(executable=lambda x: None, params=[2])
    task2: Task = Task(executable=lambda x: None, params=[3])

    task3: Task = Task(executable=lambda x: None, params=[5])
    task4: Task = Task(executable=lambda x: None, params=[7])
    task5: Task = Task(executable=lambda x: None, params=[9])
    tasklist1: TaskList = TaskList(tasks=[task3, task4, task5], chain_output=False)

    assert exe.add_task(task1) == task1.uid
    assert exe.get_task_position(task1.uid) == 0

    assert exe.add_task(task2) == task2.uid
    assert exe.get_task_position(task2.uid) == 1

    assert exe.add_task(tasklist1) == tasklist1.uid
    assert exe.get_task_position(tasklist1.uid) == 2


# test simple run one task and finish
def test_executor_run_one():
    make_directory(TEST_DIR)
    exe = ExecutorCore(terminate_on_finished=True, block_until_finished=False)

    task1: Task = Task(executable=write_file, params=["test_a", TEST_DIR, 'foo', False], name="test_executor_run_one")
    exe.add_task(task1)

    exe.run()
    while exe.is_active():
        time.sleep(0.0001)

    assert check_file('test_a', TEST_DIR, 'foo')
    delete_directory(TEST_DIR)


# test simple run one task and finish
def test_executor_run():
    make_directory(TEST_DIR)
    exe = ExecutorCore(terminate_on_finished=True, block_until_finished=False)

    task1: Task = Task(executable=write_file, params=["test_b", TEST_DIR, 'foo', False], name="test_executor_run_1a")
    exe.add_task(task1)
    task2: Task = Task(executable=write_file, params=["test_c", TEST_DIR, 'bar', False], name="test_executor_run_1b")
    exe.add_task(task2)

    def write_file_chainer(name, directory, content, append, prev) -> str:
        write_file(name, directory, prev + content, append)
        return prev + content

    task3: Task = Task(executable=write_file, params=["test_d", TEST_DIR, 'lorem', False],
                       name="test_executor_run_2a")
    task4: Task = Task(executable=write_file_chainer, params=["test_d", TEST_DIR, 'ipsum', True],
                       name="test_executor_run_2b")
    task5: Task = Task(executable=write_file_chainer, params=["test_d", TEST_DIR, 'dolor', True],
                       name="test_executor_run_2c")
    tasklist1: TaskList = TaskList(tasks=[task3, task4, task5], chain_output=True, name="test_executor_run_2")
    exe.add_task(tasklist1)

    exe.run()
    while exe.is_active():
        time.sleep(0.0001)

    time.sleep(0.1)  # time for changes to propagate
    assert check_file('test_b', TEST_DIR, 'foo')
    assert check_file('test_c', TEST_DIR, 'bar')
    assert check_file('test_d', TEST_DIR, 'loremloremipsumloremipsumdolor')

    delete_directory(TEST_DIR)


def test_executor_enqueue_during():
    make_directory(TEST_DIR)
    exe = ExecutorCore(terminate_on_finished=True, block_until_finished=False)

    task1: Task = Task(executable=write_file, params=["test_a", TEST_DIR, 'foo', False])
    exe.add_task(task1)
    task2: Task = Task(executable=write_file, params=["test_b", TEST_DIR, 'bar', False])
    exe.add_task(task2)

    exe.run()

    task3: Task = Task(executable=write_file, params=["test_c", TEST_DIR, 'boo', False])
    exe.add_task(task3)

    while exe.is_active():
        time.sleep(0.0001)

    time.sleep(0.1)  # time for changes to propagate
    assert check_file('test_a', TEST_DIR, 'foo')
    assert check_file('test_b', TEST_DIR, 'bar')
    assert check_file('test_c', TEST_DIR, 'boo')

    delete_directory(TEST_DIR)


def test_executor_fail():
    make_directory(TEST_DIR)
    exe = ExecutorCore(terminate_on_finished=True, block_until_finished=False)

    task1: Task = Task(executable=write_file, params=["test_a", TEST_DIR, 'foo', False])
    exe.add_task(task1)

    def kill():
        assert False
    task2: Task = Task(executable=kill, params=[])
    exe.add_task(task2)

    task3: Task = Task(executable=write_file, params=["test_c", TEST_DIR, 'boo', False])
    exe.add_task(task3)
    exe.run()

    while exe.is_active():
        time.sleep(0.0001)

    time.sleep(0.1)  # time for changes to propagate
    assert check_file('test_a', TEST_DIR, 'foo')
    assert check_file('test_c', TEST_DIR, 'boo')

# test singleton is working
# it is, duh