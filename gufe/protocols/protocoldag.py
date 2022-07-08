# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from typing import Iterable, List, Dict, Set, Optional

import networkx as nx

from .protocolunit import ProtocolUnit, ProtocolUnitToken
from .results import ProtocolDAGResult


class ProtocolDAG:
    """An executable directed, acyclic graph (DAG) composed of `ProtocolUnit`s
    with dependencies specified.

    This is the unit of execution passed to an alchemical `Scheduler`.
    A `ProtocolDAG` yields a `ProtocolResult`, which can be placed in a `ResultStore`.

    """

    def __init__(
        self,
        *,
        protocol_units: Iterable[ProtocolUnit],
        name: Optional[str] = None,
    ):
        """Create a new `ProtocolDAG`.

        Parameters
        ----------
        protocol_units : nx.DiGraph
            The `ProtocolUnit`s that make up this `ProtocolDAG`, with
            dependencies set as inputs.
        name : str
            Unique identifier for this `ProtocolDAG`.

        """
        self._name = name

        # build graph from protocol units
        self._graph = self._build_graph(protocol_units)

    def _build_graph(self, protocol_units):
        G = nx.DiGraph()
        """Build dependency DAG of ProtocolUnits with input keys stored on edges"""

        # build mapping of tokens to ProtocolUnits
        mapping = {pu.token: pu for pu in protocol_units}

        for pu in protocol_units:
            dependencies = []
            for key, value in pu.inputs.items():
                if isinstance(value, dict):
                    for k, v in value.items():
                        if isinstance(v, ProtocolUnitToken):
                            G.add_edge(pu, mapping[v], inputs={key: {k: None}})
                elif isinstance(value, list):
                    for i in value:
                        if isinstance(v, ProtocolUnitToken):
                            G.add_edge(pu, mapping[v], inputs={key: [None]})
                elif isinstance(value, ProtocolUnitToken):
                    G.add_edge(pu, mapping[v], inputs={key: None})

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

        # build mapping of tokens to ProtocolUnits
        mapping = {pu.token: pu for pu in graph.nodes}

        # iterate in DAG order
        for unit in reversed(list(nx.topological_sort(self._graph))):

            inputs = unit.inputs

            # for each successor, get result and input data;
            for d in self._graph.successors(unit):
                edge = self._graph.edges[unit, d]['inputs']
                result = graph.nodes[d]['result']

            ninputs = dict()
            for key, value in inputs.items():
                if isinstance(value, dict):
                    if key not in ninputs:
                        ninputs[key] = dict()
                    for k, v in value.items():
                        if isinstance(v, ProtocolUnitToken):
                            ninputs[key][k] = mapping[v]['result']
                            G.add_edge(pu, mapping[v], inputs={key: {k: None}})
                elif isinstance(value, list):
                    for i in value:
                        if isinstance(v, ProtocolUnitToken):
                            G.add_edge(pu, mapping[v], inputs={key: [None]})
                elif isinstance(value, ProtocolUnitToken):
                    G.add_edge(pu, mapping[v], inputs={key: None})
                
            # construct input to unit
            # execute
            result = unit.execute(
                inputs=[
                    graph.nodes[d]["result"]

                    for d in self._graph.successors(unit)
                ]
            )
            # attach results
            graph.nodes[unit]["result"] = result

        return ProtocolDAGResult(name=self._name, graph=graph)
