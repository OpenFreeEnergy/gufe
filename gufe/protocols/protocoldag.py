# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from typing import Iterable, List, Dict, Set, Optional, Union

import networkx as nx

from .protocolunit import ProtocolUnit, ProtocolUnitKey
from .results import ProtocolUnitResult, ProtocolDAGResult, ProtocolDAGFailure


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

        # build mapping of keys to ProtocolUnits
        mapping = {pu.key: pu for pu in protocol_units}

        for pu in protocol_units:
            for key, value in pu.inputs.items():
                if isinstance(value, dict):
                    for k, v in value.items():
                        if isinstance(v, ProtocolUnit):
                            G.add_edge(pu.to_keyed_dict(), v.to_keyed_dict())
                elif isinstance(value, list):
                    for i in value:
                        if isinstance(i, ProtocolUnit):
                            G.add_edge(pu.to_keyed_dict(), i.to_keyed_dict())
                elif isinstance(value, ProtocolUnit):
                    G.add_edge(pu.to_keyed_dict(), value.to_keyed_dict())

        return G

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
        for unit_dict in reversed(list(nx.topological_sort(self._graph))):

            unit: ProtocolUnit = ProtocolUnit.from_keyed_dict(unit_dict)

            # translate each `ProtocolUnitKey` in input into corresponding
            # `ProtocolUnitResult`
            inputs = self._pu_to_pur(unit.inputs, graph)
            
            # execute
            result = unit.execute(**inputs)

            # attach result to this `ProtocolUnit`
            graph.nodes[unit]["result"] = result.to_keyed_dict()

            if not result.ok():
                return ProtocolDAGFailure(name=self._name, graph=graph)

        return ProtocolDAGResult(name=self._name, graph=graph)

    def _pu_to_pur(
            self,
            inputs, 
            graph):
        """Convert each `ProtocolUnit` to its corresponding `ProtocolUnitResult`.

        """
        if isinstance(inputs, dict):
            return {key: self._pu_to_pur(value, graph) for key, value in inputs.items()}
        elif isinstance(inputs, list):
            return [self._pu_to_pur(value, graph) for value in inputs]
        elif isinstance(inputs, ProtocolUnit):
            return graph.nodes[inputs]['result']
        else:
            return inputs


