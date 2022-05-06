import abc
from typing import Iterable, List

from .results import ProtocolUnitResult, ProtocolResult


class ProtocolUnit(abc.ABC):
    """A unit of work computable by

    """

    def __init__(
            self,
            settings: "ProtocolSettings",
            *dependencies: Iterable["ProtocolUnit"],
        ):

        self._settings = settings
        self._dependencies = dependencies
        self._status = "WAITING"


    def execute(self, block=True):
        """Given `ProtocolUnitResult`s from dependencies, execute this `ProtocolUnit`.

        Parameters
        ----------
        block : bool
            If `True`, block until execution completes; otherwise run in its own thread.

        """
        if block:
            dep_results = [dep.estimate() for dep in self._dependencies]

            self._status = "RUNNING"
            self._execute(dep_results)
            self._status = "COMPLETE"

        else:
            #TODO: wrap in a thread; update status
            ...



    @abc.abstractmethod
    def _execute(self, dependency_results: List[ProtocolUnitResult]):
        ...

    @abc.abstractmethod
    def estimate(self) -> ProtocolUnitResult:
        ...

    def status(self):
        ...

    def get_artifacts(self) -> bytes:
        ...


class ProtocolDAG:
    """An executable directed, acyclic graph (DAG) composed of `ProtocolUnit`s
    with dependencies specified.

    This is the unit of execution passed to an alchemical `Scheduler`.
    A `ProtocolDAG` yields a `ProtocolResult`, which can be placed in a `ResultStore`.

    """

    def __init__(
            self,
            protocol_units: Iterable[ProtocolUnit]
                ):
            ...

    def to_dask(self):
        """Produce a `dask`-executable DAG from this `ProtocolDAG` as a `dask.Delayed` object.

        """
        ...

    def execute(self) -> ProtocolResult:
        """Execute the full DAG in-serial, in process.

        """
        ...
