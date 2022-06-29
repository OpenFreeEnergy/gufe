# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from typing import Iterable, List, Dict, Set, Optional, Union

import networkx as nx

from .protocolunit import ProtocolUnit
from .results import ProtocolDAGResult, ProtocolDAGFailure


class ProtocolDAG:
    """An executable directed, acyclic graph (DAG) composed of `ProtocolUnit`s
    with dependencies specified.

    This is the unit of execution passed to an alchemical `Scheduler`.
    A `ProtocolDAG` yields a `ProtocolResult`, which can be placed in a `ResultStore`.

    """

    def __init__(
        self,
        graph: nx.DiGraph,
        name: Optional[str] = None,
    ):
        """Create a new `ProtocolDAG`.

        Parameters
        ----------
        graph : nx.DiGraph
            The `ProtocolUnit`s, with dependencies set, as a networkx `DiGraph`.
        name : str
            Unique identifier for this `ProtocolDAG`.

        """
        self._graph = graph
        self._name = name

    @property
    def name(self):
        return self._name

    @property
    def graph(self):
        return self._graph

    def execute(self) -> Union[ProtocolDAGResult, ProtocolDAGFailure]:
        """Execute the full DAG in-serial, in process."""
        # operate on a copy, since we'll add ProtocolUnitResults as node attributes
        graph = self._graph.copy(as_view=False)

        # iterate in DAG order
        for unit in reversed(list(nx.topological_sort(self._graph))):
            result = unit.execute(
                dependency_results=[
                    graph.nodes[d]["result"]
                    for d in self._graph.successors(unit)
                ]
            )
            # attach results
            graph.nodes[unit]["result"] = result

            if not result.ok():
                return ProtocolDAGFailure(name=self._name, graph=graph)

        return ProtocolDAGResult(name=self._name, graph=graph)
