# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from typing import Iterable, List, Dict, Set, Optional

import networkx as nx

from .protocolunit import ProtocolUnit, ProtocolUnitToken
from .results import ProtocolUnitResult, ProtocolDAGResult


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
        protocol_units : Iterable[ProtocolUnit]
            The `ProtocolUnit`s that make up this `ProtocolDAG`, with
            dependencies included as inputs.
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
            for key, value in pu.inputs.items():
                if isinstance(value, dict):
                    for k, v in value.items():
                        if isinstance(v, ProtocolUnitToken):
                            G.add_edge(pu, mapping[v])
                elif isinstance(value, list):
                    for i in value:
                        if isinstance(i, ProtocolUnitToken):
                            G.add_edge(pu, mapping[i])
                elif isinstance(value, ProtocolUnitToken):
                    G.add_edge(pu, mapping[value])

        return G

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
        mapping = {pu.token: pu for pu in graph}

        # iterate in DAG order
        for unit in reversed(list(nx.topological_sort(self._graph))):

            # translate each `ProtocolUnitToken` in input into corresponding
            # `ProtocolUnitResult`
            inputs = self._detokenize_dependencies(unit.inputs, graph, mapping)
                
            # execute
            result = unit.execute(**inputs)

            # attach result to this `ProtocolUnit`
            graph.nodes[unit]["result"] = result

        return ProtocolDAGResult(name=self._name, graph=graph)

    @staticmethod
    def _detokenize_dependencies(
            inputs, 
            graph,
            mapping: Dict[str, ProtocolUnit]):
        ninputs = dict()
        for key, value in inputs.items():
            if isinstance(value, dict):
                if key not in ninputs:
                    ninputs[key] = dict()
                for k, v in value.items():
                    if isinstance(v, ProtocolUnitToken):
                        ninputs[key][k] = graph.nodes[mapping[v]]['result']
                    else:
                        ninputs[key][k] = v
            elif isinstance(value, list):
                if key not in ninputs:
                    ninputs[key] = list()
                for i in value:
                    if isinstance(i, ProtocolUnitToken):
                        ninputs[key].append(graph.nodes[mapping[i]]['result'])
                    else:
                        ninputs[key].append(i)
            elif isinstance(value, ProtocolUnitToken):
                ninputs[key] = graph.nodes[mapping[value]]['result']
            else:
                ninputs[key] = value

        return ninputs

