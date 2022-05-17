# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from typing import Iterable, List, Dict

from .protocolunit import ProtocolUnit
from .results import ProtocolDAGResult


class ProtocolDAG:
    """An executable directed, acyclic graph (DAG) composed of `ProtocolUnit`s
    with dependencies specified.

    This is the unit of execution passed to an alchemical `Scheduler`.
    A `ProtocolDAG` yields a `ProtocolResult`, which can be placed in a `ResultStore`.

    """

    def __init__(
            self,
            name: str,
            protocol_units: Iterable[ProtocolUnit]
        ):
        """Create a new `ProtocolDAG`.

        Parameters
        ----------
        name : str
            Name of the `Protocol` that generated this `ProtocolDAG`.
        protocol_units : Iterable[ProtocolUnit]
            The `ProtocolUnit`s, with dependencies set, to include in the `ProtocolDAG`.

        """
        self._name = name
        self._protocol_units = tuple(protocol_units)

    @property
    def name(self):
        return self._name

    def execute(self) -> ProtocolDAGResult:
        """Execute the full DAG in-serial, in process.

        """
        completed: List[ProtocolUnit] = []
        while len(completed) != len(self._protocol_units):
            for pu in self._protocol_units:
                if pu.status == 'READY':
                    pu.execute()
                    completed.append(pu)

        return ProtocolDAGResult(name=self._name, units=[pu.result for pu in completed])
