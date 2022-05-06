import abc
from typing import Iterable, Union

from .scheduler import Scheduler
from .results import ResultStore

class Executor(abc.ABC):
    ...
    def __init__(self,
            schedulers: Union[Scheduler, Iterable[Scheduler]],
            results: ResultStore):
        ...


class LocalExecutor(Executor):
    """In-process, synchronous implementation of an `Executor`.

    Useful for testing and single-machine execution.

    """
    ...
