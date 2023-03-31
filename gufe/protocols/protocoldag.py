# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from copy import copy
from collections import defaultdict
import os
from typing import Iterable, Optional, Union, Any
from os import PathLike
from pathlib import Path
import tempfile

import networkx as nx

from ..tokenization import GufeTokenizable, GufeKey
from .protocolunit import (
    ProtocolUnit, ProtocolUnitResult, ProtocolUnitFailure, Context
)


class DAGMixin:
    _protocol_units: list[ProtocolUnit]

    _name: Optional[str]
    _graph: nx.DiGraph

    # labels for identifying source of this DAG

    ## key of the Transformation that this DAG corresponds to
    _transformation_key: Union[GufeKey, None]

    ## key of the ProtocolDAG this DAG extends
    _extends_key: Optional[GufeKey]

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
        """DAG of `ProtocolUnit`s that comprise this object.

        """
        return self._graph

    @property
    def protocol_units(self) -> list[ProtocolUnit]:
        """List of `ProtocolUnit`s given in DAG-order.

        """
        return list(self._iterate_dag_order(self._graph))

    @property
    def transformation_key(self) -> Union[GufeKey, None]:
        """The `GufeKey` of the `Transformation` this object performs.

        If `None`, then this object was not created from a `Transformation`.
        This may be the case when creating a `ProtocolDAG` from a `Protocol`
        directly, without use of a `Transformation` object.

        This functions as a label, indicating where this object came from.

        """
        return self._transformation_key

    @property
    def extends_key(self) -> Union[GufeKey, None]:
        """The `GufeKey` of the `ProtocolDAGResult` this object extends.

        If `None`, then this object does not extend from a result at all.

        This functions as a label, indicating where this object came from.
        It can be used to reconstruct the set of extension relationships
        between a collection of ProtocolDAGs.

        """
        return self._extends_key


class ProtocolDAGResult(GufeTokenizable, DAGMixin):
    """Result for a single `ProtocolDAG` execution.

    There may be many of these for a given `Transformation`. Data elements from
    these objects are combined by `Protocol.gather` into a `ProtocolResult`.

    Attributes
    ----------
    name : str
        Optional identifier for this `ProtocolDAGResult`.
    protocol_units : list[ProtocolUnit]
        `ProtocolUnit`s (given in DAG-dependency order) used to compute this
        `ProtocolDAGResult`.
    protocol_unit_results : list[ProtocolUnitResult]
        `ProtocolUnitResult`s (given in DAG-dependency order) corresponding to
        each `ProtocolUnit` used to compute this `ProtocolDAGResult`.
    graph : nx.DiGraph
        Graph of `ProtocolUnit`s as nodes, with directed edges to each
        `ProtocolUnit`'s dependencies.
    result_graph : nx.DiGraph
        Graph of `ProtocolUnitResult`s as nodes, with directed edges to each
        `ProtocolUnitResult`'s dependencies.
    transformation_key : Union[GufeKey, None]
        Key of the `Transformation` that this `ProtocolDAGResult` corresponds
        to, if applicable. This functions as a label for identifying the source
        of this `ProtocolDAGResult`.
    extends_key : Optional[GufeKey]
        Key of the `ProtocolDAGResult` that this `ProtocolDAGResult` extends from.
        This functions as a label for identifying the source of this `ProtocolDAGResult`;
        it can be used to reconstruct the tree of extensions from a collection
        of `ProtocolUnitResult`\s.

    """
    _protocol_unit_results: list[ProtocolUnitResult]
    _unit_result_mapping: dict[ProtocolUnit, list[ProtocolUnitResult]]
    _result_unit_mapping: dict[ProtocolUnitResult, ProtocolUnit]


    def __init__(
        self, 
        *,
        protocol_units: list[ProtocolUnit],
        protocol_unit_results: list[ProtocolUnitResult],
        transformation_key: Union[GufeKey, None],
        extends_key: Optional[GufeKey] = None,
        name: Optional[str] = None,
    ):
        self._name = name
        self._protocol_units = protocol_units
        self._protocol_unit_results = protocol_unit_results

        self._transformation_key = GufeKey(transformation_key) if transformation_key is not None else None
        self._extends_key = GufeKey(extends_key) if extends_key is not None else None

        # build graph from protocol units
        self._graph = self._build_graph(protocol_units)

        # build graph from protocol unit results
        self._result_graph = self._build_graph(protocol_unit_results)

        # build mapping from protocol units to results
        keys_to_pu = {unit.key: unit for unit in self._protocol_units}
        unit_result_mapping = defaultdict(list)
        self._result_unit_mapping = dict()
        for result in protocol_unit_results:
            pu = keys_to_pu[result.source_key]
            unit_result_mapping[pu].append(result)
            self._result_unit_mapping[result] = pu
        self._unit_result_mapping = dict(unit_result_mapping)

    @classmethod
    def _defaults(cls):
        # not used by `ProtocolDAG`
        return {}

    def _to_dict(self):
        return {'name': self.name,
                'protocol_units': self._protocol_units,
                'protocol_unit_results': self._protocol_unit_results,
                'transformation_key': self._transformation_key,
                'extends_key': self._extends_key}

    @classmethod
    def _from_dict(cls, dct: dict):
        return cls(**dct)

    @property
    def result_graph(self) -> nx.DiGraph:
        return self._result_graph

    @property
    def protocol_unit_results(self) -> list[ProtocolUnitResult]:
        return list(self._iterate_dag_order(self.result_graph))

    @property
    def protocol_unit_failures(self) -> list[ProtocolUnitFailure]:
        """A list of all failed units.

        Note
        ----
        These are returned in DAG order
        """
        # mypy can't figure out the types here, .ok() will ensure a certain type
        # https://mypy.readthedocs.io/en/stable/common_issues.html?highlight=cast#complex-type-tests
        return [r for r in self.protocol_unit_results if not r.ok()]  # type: ignore
    
    @property
    def protocol_unit_successes(self) -> list[ProtocolUnitResult]:
        """A list of only successful `ProtocolUnit` results.

        Note
        ----
        These are returned in DAG order
        """
        return [r for r in self.protocol_unit_results if r.ok()]

    def unit_to_result(self, protocol_unit: ProtocolUnit) -> ProtocolUnitResult:
        """Return the successful result for a given Unit.

        Returns
        -------
        success : ProtocolUnitResult
          the successful result for this Unit

        Raises
        ------
        KeyError
          if either there are no results, or only failures
        """
        try:
            units = self._unit_result_mapping[protocol_unit]
        except KeyError:
            raise KeyError("No such `protocol_unit` present")
        else:
            for u in units:
                if u.ok():
                    return u
            else:
                raise KeyError("No success for `protocol_unit` found")

    def unit_to_all_results(self, protocol_unit: ProtocolUnit) -> list[ProtocolUnitResult]:
        """Return all results (sucess and failure) for a given Unit.

        Returns
        -------
        results : list[ProtocolUnitResult]
          results for a given unit

        Raises
        ------
        KeyError
          if no results present for a given unit
        """
        try:
            return self._unit_result_mapping[protocol_unit]
        except KeyError:
            raise KeyError("No such `protocol_unit` present")

    def result_to_unit(self, protocol_unit_result: ProtocolUnitResult) -> ProtocolUnit:
        try:
            return self._result_unit_mapping[protocol_unit_result]
        except KeyError:
            raise KeyError("No such `protocol_unit_result` present")

    def ok(self) -> bool:
        # ensure that for every protocol unit, there is an OK result object
        return all(any(pur.ok() for pur in self._unit_result_mapping[pu])
                   for pu in self._protocol_units)

    @property
    def terminal_protocol_unit_results(self) -> list[ProtocolUnitResult]:
        """Get ProtocolUnitResults that terminate the DAG.

        Returns
        -------
        list[ProtocolUnitResult]
          All ProtocolUnitResults which do not have a ProtocolUnitResult that
          follows on (depends) on them.
        """
        return [u for u in self._protocol_unit_results
                if not nx.ancestors(self._result_graph, u)]


class ProtocolDAG(GufeTokenizable, DAGMixin):
    """An executable directed, acyclic graph (DAG) composed of `ProtocolUnit`
    with dependencies specified.

    A single `ProtocolDAG` execution should yield sufficient information to
    calculate a free energy difference (though perhaps not converged) between
    two `ChemicalSystem` objects.
    
    A `ProtocolDAG` yields a `ProtocolDAGResult` when executed.

    Attributes
    ----------
    name : str
        Optional identifier for this `ProtocolDAGResult`.
    protocol_units : list[ProtocolUnit]
        `ProtocolUnit`s (given in DAG-dependency order) used to compute this
        `ProtocolDAGResult`.
    graph : nx.DiGraph
        Graph of `ProtocolUnit`s as nodes, with directed edges to each
        `ProtocolUnit`'s dependencies.
    transformation_key : Union[GufeKey, None]
        Key of the `Transformation` that this `ProtocolDAG` corresponds to, if
        applicable. This functions as a label for identifying the source of
        this `ProtocolDAG`. This label will be passed on to the
        `ProtocolDAGResult` resulting from execution of this `ProtocolDAG`.
    extends_key : Optional[GufeKey]
        Key of the `ProtocolDAGResult` that this `ProtocolDAG` extends from.
        This functions as a label for identifying the source of this
        `ProtocolDAG`. This label will be passed on to the
        `ProtocolDAGResult` resulting from execution of this `ProtocolDAG`.

    """

    def __init__(
        self,
        *,
        protocol_units: list[ProtocolUnit],
        transformation_key: Union[GufeKey, None],
        extends_key: Optional[GufeKey] = None,
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

        self._transformation_key = GufeKey(transformation_key) if transformation_key is not None else None
        self._extends_key = GufeKey(extends_key) if extends_key is not None else None

        # build graph from protocol units
        self._graph = self._build_graph(protocol_units)

    @classmethod
    def _defaults(cls):
        # not used by `ProtocolDAG`
        return {}

    def _to_dict(self):
        return {'name': self.name,
                'protocol_units': self.protocol_units,
                'transformation_key': self._transformation_key,
                'extends_key': self._extends_key}

    @classmethod
    def _from_dict(cls, dct: dict):
        return cls(**dct)


def execute_DAG(protocoldag: ProtocolDAG, *,
                shared: Path,
                scratch_basedir: Path,
                keep_scratch: bool = False,
                raise_error: bool = True,
                ) -> ProtocolDAGResult:
    """Execute the full DAG in-serial, in process.

    This is intended for debug use for Protocol developers.
    Running locally is generally not useful for production purposes.

    Parameters
    ----------
    protocoldag : ProtocolDAG
        The `ProtocolDAG` to execute.
    shared : Path
        Path to scratch space that persists across whole DAG execution, but is
        removed after. Used by some `ProtocolUnit`s to pass file contents to
        dependent `ProtocolUnit`s. If not given, defaults to os cwd (current
        directory).
    scratch_basedir : Path
        Filesystem path to use for `ProtocolUnit` `scratch` space.
    keep_scratch : bool
        If True, don't remove scratch directories for `ProtocolUnit`s after
        completion.
    raise_error : bool
        If True, raise an exception if a ProtocolUnit fails, default True
        if False, any exceptions will be stored as `ProtocolUnitFailure`
        objects inside the returned `ProtocolDAGResult`

    Returns
    -------
    ProtocolDAGResult
        The result of executing the `ProtocolDAG`.

    """
    # iterate in DAG order
    results: dict[GufeKey, ProtocolUnitResult] = {}
    for unit in protocoldag.protocol_units:

        # translate each `ProtocolUnit` in input into corresponding
        # `ProtocolUnitResult`
        inputs = _pu_to_pur(unit.inputs, results)

        scratch_tmp = tempfile.TemporaryDirectory(
                prefix=f"{str(unit.key)}__",
                dir=scratch_basedir)
        scratch = Path(scratch_tmp.name)

        context = Context(shared=shared,
                          scratch=scratch)

        # execute
        result = unit.execute(
                context=context,
                raise_error=raise_error,
                **inputs)

        if not keep_scratch:
            scratch_tmp.cleanup()

        # attach result to this `ProtocolUnit`
        results[unit.key] = result

        if not result.ok():
            break

    return ProtocolDAGResult(
            name=protocoldag.name, 
            protocol_units=protocoldag.protocol_units, 
            protocol_unit_results=list(results.values()),
            transformation_key=protocoldag.transformation_key,
            extends_key=protocoldag.extends_key)


def _pu_to_pur(
        inputs: Union[dict[str, Any], list[Any], ProtocolUnit],
        mapping: dict[GufeKey, ProtocolUnitResult]):
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

