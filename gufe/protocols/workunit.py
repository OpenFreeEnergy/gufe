import abc
from typing import List

class WorkUnit(abc.ABC):
    """A unit of work, gives status of units.

    """

    def __init__(
            self,
            ):
        ...

    def execute(self):
        ...

    def estimate(self) -> Result:
        ...

    def status(self):
        ...

    def get_artifacts(self) -> bytes:
        ...
