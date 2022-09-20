# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from typing import Iterable, List, Dict, Set, Optional, Union, Any
from os import PathLike
from pathlib import Path
import tempfile

import networkx as nx

from ..tokenization import GufeTokenizable, GufeKey
from .protocolunit import ProtocolUnit, ProtocolUnitResult, ProtocolUnitResultBase


class DAGMixin:
    _name: Optional[str]
    _graph: nx.DiGraph

    @staticmethod 
    def _build_graph(nodes):
        """Build dependency DAG of ProtocolUnits with input keys stored on edges"""
        G = nx.DiGraph()

        for node in nodes:
            G.add_node(node)
            for dep in node.dependencies:
                G.add_edge(node, dep)

        return G

    @staticmethod
    def _iterate_dag_order(graph):
        return reversed(list(nx.topological_sort(graph)))

    @property
    def name(self):
        return self._name

    @property
    def graph(self) -> nx.DiGraph:
        """DAG of `ProtocolUnit`s that produced this `ProtocolDAGResult`.

        """
        return self._graph
    @property
    def protocol_units(self):
        """List of `ProtocolUnit`s given in DAG-order.

        """
        return list(self._iterate_dag_order(self._graph))


class ProtocolDAGResult(GufeTokenizable, DAGMixin):
    """Result for a single `ProtocolDAG` execution.

    There may be many of these for a given `Transformation`. Data elements from
    these objects are combined by `Protocol.gather` into a `ProtocolResult`.

    Attributes
    ----------
    name : str
        Optional identifier for this `ProtocolDAGResult`.
    protocol_units : List[ProtocolUnit]
        `ProtocolUnit`s (given in DAG-dependency order) used to compute this
        `ProtocolDAGResult`.
    protocol_unit_results : List[ProtocolUnit]
        `ProtocolUnitResult`s (given in DAG-dependency order) corresponding to
        each `ProtocolUnit` used to compute this `ProtocolDAGResult`.
    graph : nx.DiGraph
        Graph of `ProtocolUnit`s as nodes, with directed edges to each
        `ProtocolUnit`'s dependencies.
    result_graph : nx.DiGraph
        Graph of `ProtocolUnitResult`s as nodes, with directed edges to each
        `ProtocolUnitResult`'s dependencies.

    """
    def __init__(self, *, 
            name=None, 
            protocol_units: List[ProtocolUnit], 
            protocol_unit_results: List[ProtocolUnitResult]):
        self._name = name
        self._protocol_units = protocol_units
        self._protocol_unit_results = protocol_unit_results

        # build graph from protocol units
        self._graph = self._build_graph(protocol_units)

        # build graph from protocol unit results
        self._result_graph = self._build_graph(protocol_unit_results)

        # build mapping from protocol units to results
        # TODO: currently assumes only a single ProtocolUnitResult for each
        # ProtocolUnit; will need to update this class if multiple failed
        # ProtocolUnitResults for a given ProtocolUnitResult supported in this
        # data structure
        keys_to_pu = {unit.key: unit for unit in self._protocol_units}
        self._unit_result_mapping = {keys_to_pu[result.source_key]: result
                                     for result in self._protocol_unit_results}
        self._result_unit_mapping = {result: unit for unit, result
                                     in self._unit_result_mapping.items()}

    def _defaults(self):
        # not used by `ProtocolDAG`
        return {}

    def _to_dict(self):
        return {'name': self.name,
                'protocol_units': self._protocol_units,
                'protocol_unit_results': self._protocol_unit_results}

    @classmethod
    def _from_dict(cls, dct: Dict):
        return cls(**dct)

    @property
    def result_graph(self) -> nx.DiGraph:
        """DAG of `ProtocolUnitResult`s that compose this `ProtocolDAGResult`.

        """
        return self._result_graph

    @property
    def protocol_unit_results(self):
        return list(self._iterate_dag_order(self.result_graph))

    @property
    def protocol_unit_failures(self):
        """A list of `ProtocolUnitFailure`s corresponding to only failed
        `ProtocolUnit`s.

        """
        return [r for r in self.protocol_unit_results if not r.ok()]

    def unit_to_result(self, protocol_unit: ProtocolUnit):
        try:
            return self._unit_result_mapping[protocol_unit]
        except KeyError:
            raise KeyError("No such `protocol_unit` present")

    def result_to_unit(self, protocol_unit_result: ProtocolUnitResult):
        try:
            return self._result_unit_mapping[protocol_unit_result]
        except KeyError:
            raise KeyError("No such `protocol_unit_result` present")

    def ok(self) -> bool:
        # ensure that for every protocol unit, there is definitely a
        # corresponding result
        return set(self._protocol_units) == set(self._unit_result_mapping.keys())


class ProtocolDAG(GufeTokenizable, DAGMixin):
    """An executable directed, acyclic graph (DAG) composed of `ProtocolUnit`s
    with dependencies specified.

    A single `ProtocolDAG` execution should yield sufficient information to
    calculate a free energy difference (though perhaps not converged) between
    two `ChemicalSystem`s.
    
    A `ProtocolDAG` yields a `ProtocolDAGResult` when executed.

    Attributes
    ----------
    name : str
        Optional identifier for this `ProtocolDAGResult`.
    protocol_units : List[ProtocolUnit]
        `ProtocolUnit`s (given in DAG-dependency order) used to compute this
        `ProtocolDAGResult`.
    graph : nx.DiGraph
        Graph of `ProtocolUnit`s as nodes, with directed edges to each
        `ProtocolUnit`'s dependencies.

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


def execute(protocoldag: ProtocolDAG, *, 
            shared: PathLike = None) -> ProtocolDAGResult:
    """Execute the full DAG in-serial, in process.

    This is intended for debug use for Protocol developers.
    Running locally is generally not useful for production purposes.

    Parameters
    ----------
    protocoldag : ProtocolDAG
        The `ProtocolDAG` to execute.
    shared : Optional[PathLike]
       Path to scratch space that persists across whole DAG execution, but
       is removed after. Used by some `ProtocolUnit`s to pass file contents
       to dependent `ProtocolUnit`s.

    Returns
    -------
    ProtocolDAGResult
        The result of executing the `ProtocolDAG`.

    """
    if shared is None:
        shared_tmp = tempfile.TemporaryDirectory()
        shared_ = Path(shared_tmp.name)
    else:
        shared_ = Path(shared)

    # iterate in DAG order
    results: Dict[GufeKey, ProtocolUnitResult] = {}
    for unit in protocoldag.protocol_units:

        # translate each `ProtocolUnit` in input into corresponding
        # `ProtocolUnitResult`
        inputs = _pu_to_pur(unit.inputs, results)

        # execute
        result = unit.execute(shared=shared_, **inputs)

        # attach result to this `ProtocolUnit`
        results[unit.key] = result

        if not result.ok():
            break

    # TODO: change this part once we have clearer ideas on how to inject
    # persistent storage use
    if shared is None:
        shared_tmp.cleanup()

    return ProtocolDAGResult(
            name=protocoldag.name, 
            protocol_units=protocoldag.protocol_units, 
            protocol_unit_results=list(results.values()))


def _pu_to_pur(
        inputs: Union[Dict[str, Any], List[Any], ProtocolUnit],
        mapping: Dict[GufeKey, ProtocolUnitResult]):
    """Convert each `ProtocolUnit` found within `inputs` to its corresponding
    `ProtocolUnitResult`.

    Parameters
    ----------
    inputs
        Arbitrarily-nested dict or list, with `ProtocolUnit`s present among
        values/elements. Can also be a single `ProtocolUnit`.

    Returns
    -------
    Data structure identical to `inputs`, except with each `ProtocolUnit`
    replaced with its corresponding `ProtocolUnitResult`.

    """
    if isinstance(inputs, dict):
        return {key: _pu_to_pur(value, mapping) for key, value in inputs.items()}
    elif isinstance(inputs, list):
        return [_pu_to_pur(value, mapping) for value in inputs]
    elif isinstance(inputs, ProtocolUnit):
        return mapping[inputs.key]
    else:
        return inputs

