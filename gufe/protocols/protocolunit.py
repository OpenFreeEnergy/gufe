# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

"""The `ProtocolUnit` class should be subclassed for all units to be used as
part of a `ProtocolDAG`.

"""

import abc
from typing import Iterable, List, Dict, Any
from pathlib import Path

from .results import ProtocolUnitResult


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
            dep_results = [dep.result for dep in self._dependencies]

            self._status = "RUNNING"
            out = self._execute(dep_results)
            self._status = "COMPLETE"
            self._result = ProtocolUnitResult(
                                name=self.__class__.__name__,
                                dependencies=dep_results, 
                                **out)

        else:
            #TODO: wrap in a thread; update status
            ...

        return self._result


    @abc.abstractmethod
    def _execute(self, dependency_results: List[ProtocolUnitResult]) -> Dict[str, Any]:
        ...

    @property
    def result(self) -> ProtocolUnitResult:
        """Return `ProtocolUnitResult` for this `ProtocolUnit`.

        Requires `status` == "COMPLETE"; exception raised otherwise.

        """
        return self._result

    @property
    def status(self):
        if self._status == "WAITING":
            # check dependencies; if all COMPLETE, change to READY
            if set(dep.status for dep in self._dependencies) == {"COMPLETE"}:
                self._status = 'READY'

        return self._status

    def get_artifacts(self) -> Dict[str, Path]:
        ...


