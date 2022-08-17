# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from typing import Iterable, List, Dict, Set, Optional, Union
from os import PathLike
import tempfile

import networkx as nx

from ..tokenize import GufeTokenizable
from .protocolunit import ProtocolUnit, ProtocolUnitResult


class ProtocolDAGResult(GufeTokenizable):
    """Result for a single `ProtocolDAG` execution.

    There may be many of these in a given `ResultStore` for a given `Transformation`.
    Data elements from these objects are combined by `Protocol.gather` into a
    `ProtocolResult`.

    Attributes
    ----------
    name : str
        Unique identifier for this `ProtocolDAGResult`.
    graph : nx.DiGraph
        The `ProtocolUnit`s, with dependencies set, as a networkx `DiGraph`.
        Each `ProtocolUnit` features its `ProtocolUnitCompletion` as a `result` attribute.

    """
    name: Optional[str]
    graph: nx.DiGraph

    def __init__(self, *, name=None, graph):
        self._name = name
        self._graph = graph

    def _defaults(self):
        # not used by `ProtocolDAG`
        return {}

    def _to_dict(self):
        return {'name': self.name,
                'graph': self.graph}

    @classmethod
    def _from_dict(cls, dct: Dict):
        return cls(**dct)

    @property
    def name(self):
        return self._name

    @property
    def graph(self):
        return self._graph

    @property
    def protocol_units(self):
        return [pu for pu in self.graph.nodes]

    @property
    def protocol_unit_results(self):
        return list(nx.get_node_attributes(self.graph, "result").values())

    def ok(self) -> bool:
        return True


class ProtocolDAGFailure(ProtocolDAGResult):

    def ok(self) -> bool:
        return False

    @property
    def protocol_unit_failures(self):
        return [r for r in nx.get_node_attributes(self.graph, "result").values() if not r.ok()]

    def _defaults(self):
        # not used by `ProtocolDAG`
        return {}

    def _to_dict(self):
        return {'name': self.name,
                'graph': self.graph}

    @classmethod
    def _from_dict(cls, dct: Dict):
        return cls(**dct)


class ProtocolDAG(GufeTokenizable):
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
        self._protocol_units = protocol_units

        # build graph from protocol units
        self._graph = self._build_graph(protocol_units)

    def _defaults(self):
        # not used by `ProtocolDAG`
        return {}

    def _to_dict(self):
        return {'name': self.name,
                'protocol_units': self.protocol_units}

    @classmethod
    def _from_dict(cls, dct: Dict):
        return cls(**dct)

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
                            G.add_edge(pu, v)
                elif isinstance(value, list):
                    for i in value:
                        if isinstance(i, ProtocolUnit):
                            G.add_edge(pu, i)
                elif isinstance(value, ProtocolUnit):
                    G.add_edge(pu, value)

        return G

    @property
    def name(self):
        return self._name

    @property
    def graph(self):
        return self._graph

    def protocol_units(self):
        return list(self._protocol_units)

    def execute(self, *, 
            dag_scratch: PathLike = None) -> Union[ProtocolDAGResult, ProtocolDAGFailure]:
        """Execute the full DAG in-serial, in process.

        Parameters
        ----------
        dag_scratch : Optional[PathLike]
           Path to scratch space that persists across whole DAG execution, but
           is removed after. Used by some `ProtocolUnit`s to pass file contents
           to dependent `ProtocolUnit`s.

        """
        if dag_scratch is None:
            dag_scratch_tmp = tempfile.TemporaryDirectory()
            dag_scratch_ = dag_scratch_tmp.name
        else:
            dag_scratch_ = dag_scratch

        # operate on a copy, since we'll add ProtocolUnitResults as node attributes
        graph = self._graph.copy(as_view=False)

        # iterate in DAG order
        for unit in reversed(list(nx.topological_sort(self._graph))):

            #unit: ProtocolUnit = ProtocolUnit.from_keyed_dict(unit_dict)

            # translate each `ProtocolUnit` in input into corresponding
            # `ProtocolUnitResult`
            inputs = self._pu_to_pur(unit.inputs, graph)

            # execute
            result = unit.execute(dag_scratch=dag_scratch_, **inputs)

            # attach result to this `ProtocolUnit`
            graph.nodes[unit]["result"] = result

            if not result.ok():
                if dag_scratch is None:
                    dag_scratch_tmp.cleanup()

                return ProtocolDAGFailure(name=self._name, graph=graph)

        # TODO: change this part once we have clearer ideas on how to inject
        # persistent storage use
        if dag_scratch is None:
            dag_scratch_tmp.cleanup()

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
