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

        return ProtocolDAGResult(units=[pu.result for pu in completed])
