# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from typing import Iterable, List, Dict, Set, Optional

import networkx as nx

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

    def execute(self) -> ProtocolDAGResult:
        """Execute the full DAG in-serial, in process."""

        # operate on a copy, since we'll add ProtocolUnitResults as node attributes
        graph = self._graph.copy(as_view=False)

        completed: Set[ProtocolUnit] = set()

        while len(completed) != len(self._graph.nodes):
            for pu in graph.nodes:
                # skip if we already completed this one
                if pu in completed:
                    continue

                dependency_edges = graph.edges(pu)

                # if all dependencies are completed, execute
                if all(
                    (
                        dependency_edge[1] in completed
                        for dependency_edge in dependency_edges
                    )
                ):
                    pur = pu.execute(
                        dependency_results=[
                            graph.nodes[dependency_edge[1]]["result"]
                            for dependency_edge in dependency_edges
                        ]
                    )

                    # attach results, add to completed list
                    graph.nodes[pu]["result"] = pur
                    completed.add(pu)

        return ProtocolDAGResult(name=self._name, graph=graph)
