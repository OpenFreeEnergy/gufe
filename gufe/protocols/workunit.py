import abc
from typing import List

from .results import WorkUnitResult

class WorkUnit(abc.ABC):
    """A unit of work computable by

    """

    def __init__(
            self,
            ):
        ...

    def execute(self):
        ...

    def estimate(self) -> WorkUnitResult:
        ...

    def status(self):
        ...

    def get_artifacts(self) -> bytes:
        ...
