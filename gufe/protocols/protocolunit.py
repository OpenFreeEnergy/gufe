# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from pathlib import Path
from typing import Iterable, List, Dict

from .results import ProtocolUnitResult, ProtocolDAGResult


class ProtocolUnit(abc.ABC):
    """A unit of work computable by

    """

    def __init__(
            self,
            settings: "ProtocolSettings",
            *dependencies: Iterable["ProtocolUnit"],
            **kwargs
        ):

        self._settings = settings
        self._dependencies = dependencies

        if self._dependencies:
            self._status = "WAITING"
        else:
            self._status = "READY"

        self._kwargs = kwargs

    @property
    def settings(self):
        return self._settings

    def execute(self, block=True) -> ProtocolUnitResult:
        """Given `ProtocolUnitResult`s from dependencies, execute this `ProtocolUnit`.

        Parameters
        ----------
        block : bool
            If `True`, block until execution completes; otherwise run in its own thread.

        """
        if block:
            dep_results = [dep.results for dep in self._dependencies]

            self._status = "RUNNING"
            out = self._execute(dep_results)
            self._status = "COMPLETE"

        else:
            #TODO: wrap in a thread; update status
            ...

        return out


    @abc.abstractmethod
    def _execute(self, dependency_results: List[ProtocolUnitResult]) -> ProtocolUnitResult:
        ...

    @property
    def results(self) -> ProtocolUnitResult:
        """Return `ProtocolUnitResult` for this `ProtocolUnit`.

        Requires `status` == "COMPLETE"; exception raised otherwise.

        """
        return self._results()

    @abc.abstractmethod
    def _results(self) -> ProtocolUnitResult:
        ...

    @property
    def status(self):
        if self._status == "WAITING":
            # check dependencies; if all COMPLETE, change to READY
            if set(dep.status for dep in self._dependencies) == {"COMPLETE"}:
                self._status = 'READY'

        return self._status

    def get_artifacts(self) -> Dict[str, Path]:
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
            self._protocol_units = tuple(protocol_units)

    def to_dask(self):
        """Produce a `dask`-executable DAG from this `ProtocolDAG` as a `dask.Delayed` object.

        """
        ...

    def execute(self) -> ProtocolDAGResult:
        """Execute the full DAG in-serial, in process.

        """
        completed = []
        while len(completed) != len(self._protocol_units):
            for pu in self._protocol_units:
                if pu.status == 'READY':
                    pu.execute()
                    completed.append(pu)

        return ProtocolDAGResult(units=[pu.results for pu in completed])
