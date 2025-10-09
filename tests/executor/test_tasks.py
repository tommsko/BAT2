from .. import delete_file, check_file, write_file, delete_directory, make_directory, TEST_DIR


from ...src.executor.task import Task, TaskList


def test_task_success():
    make_directory(TEST_DIR)

    task: Task = Task(executable=write_file, params=["test", TEST_DIR, "foo"], name="test")
    assert task.name == "test"

    assert task.execute() == "foo"
    assert task.finished_successfully

    assert check_file("test", TEST_DIR, "foo")
    delete_directory(TEST_DIR)


def test_task_fail():

    def set_checker() -> None:
        assert False, "expected assert"

    task: Task = Task(executable=set_checker, params=[], name="test")

    assert task.execute() == None
    assert not task.finished_successfully


def test_tasklist_success_isolated():
    make_directory(TEST_DIR)

    task_1: Task = Task(executable=write_file, params=['test1', TEST_DIR, 'foo', False], name="test1")
    task_2: Task = Task(executable=write_file, params=['test2', TEST_DIR, 'bar', False], name="test2")
    task_3: Task = Task(executable=write_file, params=['test3', TEST_DIR, 'oob', False], name="test3")

    tasklist = TaskList(tasks=[task_1, task_2, task_3], chain_output=False, name="list1")
    assert tasklist.name == "list1"

    assert tasklist.execute() == "oob"  # last return
    assert tasklist.finished_successfully

    assert check_file("test1", TEST_DIR, "foo")
    assert check_file("test2", TEST_DIR, "bar")
    assert check_file("test3", TEST_DIR, "oob")

    delete_directory(TEST_DIR)


def test_tasklist_success_chain_output():

    def write_file_chainer(name, directory, content, append, prev) -> str:
        write_file(name, directory, prev + content, append)
        return prev + content

    make_directory(TEST_DIR)

    task_1: Task = Task(executable=write_file, params=['test1', TEST_DIR, 'foo', False], name="test1")
    task_2: Task = Task(executable=write_file_chainer, params=['test2', TEST_DIR, 'bar', False], name="test2")
    task_3: Task = Task(executable=write_file_chainer, params=['test3', TEST_DIR, 'oob', False], name="test3")

    tasklist = TaskList(tasks=[task_1, task_2, task_3], chain_output=True, name="list1")

    assert tasklist.execute() == "foobaroob"  # last return
    assert tasklist.finished_successfully

    assert check_file("test1", TEST_DIR, "foo")
    assert check_file("test2", TEST_DIR, "foobar")
    assert check_file("test3", TEST_DIR, "foobaroob")

    delete_directory(TEST_DIR)


def test_tasklist_success_chain_output2():

    def write_file_chainer(name, directory, content, append, prev) -> str:
        new_num = str(int(prev) + 1)
        write_file(name, directory, new_num, append)
        return new_num

    make_directory(TEST_DIR)

    task_1: Task = Task(executable=write_file, params=['test', TEST_DIR, '0', False], name="test1")

    relay_tasks = [Task(executable=write_file_chainer, params=['test', TEST_DIR, None, True], name=f"test{i}")
                   for i in range(1, 10)]

    tasklist = TaskList(tasks=[task_1] + relay_tasks, chain_output=True, name="list1")

    assert tasklist.execute() == "9"
    assert tasklist.finished_successfully

    assert check_file("test", TEST_DIR, "0123456789")
    delete_directory(TEST_DIR)


def test_tasklist_fail():
    def relay(what: int) -> int:
        return what

    def fail() -> None:
        assert False, "expected assert"

    task_1: Task = Task(executable=relay, params=[0], name="test1")
    task_2: Task = Task(executable=fail, params=[], name="test2")
    task_3: Task = Task(executable=relay, params=[2], name="test3")

    tasklist = TaskList(tasks=[task_1, task_2, task_3], chain_output=False, name="list1")

    assert tasklist.execute() == None
    assert not tasklist.finished_successfully

