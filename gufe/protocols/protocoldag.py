# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
import os
import shutil
from collections import defaultdict
from collections.abc import Iterable
from copy import copy
from os import PathLike
from pathlib import Path
from typing import Any, Optional, Union

import networkx as nx

from ..tokenization import GufeKey, GufeTokenizable
from .errors import MissingUnitResultError, ProtocolUnitFailureError
from .protocolunit import Context, ProtocolUnit, ProtocolUnitFailure, ProtocolUnitResult


class DAGMixin:
    _protocol_units: list[ProtocolUnit]

    _name: str | None
    _graph: nx.DiGraph

    # labels for identifying source of this DAG

    ## key of the Transformation that this DAG corresponds to
    _transformation_key: GufeKey | None

    ## key of the ProtocolDAG this DAG extends
    _extends_key: GufeKey | None

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
        return reversed(list(nx.lexicographical_topological_sort(graph, key=lambda pu: pu.key)))

    @property
    def name(self) -> str | None:
        """Optional identifier."""
        return self._name

    @property
    def graph(self) -> nx.DiGraph:
        """DAG of ``ProtocolUnit`` nodes with edges denoting dependencies."""
        return self._graph

    @property
    def protocol_units(self) -> list[ProtocolUnit]:
        """
        List of `ProtocolUnit` s given in DAG-dependency order.

        DAG-dependency order guarantees that any task is listed after all of its
        dependencies.
        """
        return list(self._iterate_dag_order(self._graph))

    @property
    def transformation_key(self) -> GufeKey | None:
        """
        The `GufeKey` of the `Transformation` this object performs.

        If `None`, then this object was not created from a `Transformation`.
        This may be the case when creating a `ProtocolDAG` from a `Protocol`
        directly, without use of a `Transformation` object.

        This functions as a label, indicating where this object came from.
        """
        return self._transformation_key

    @property
    def extends_key(self) -> GufeKey | None:
        """The `GufeKey` of the `ProtocolDAGResult` this object extends.

        If `None`, then this object does not extend from a result at all.

        This functions as a label, indicating where this object came from.
        It can be used to reconstruct the set of extension relationships
        between a collection of ProtocolDAGs.

        """
        return self._extends_key


class ProtocolDAGResult(GufeTokenizable, DAGMixin):
    """
    Result for a single execution of an entire :class:`ProtocolDAG`.

    There may be many of these for a given `Transformation`. Data elements from
    these objects are combined by `Protocol.gather` into a `ProtocolResult`.
    """

    _protocol_unit_results: list[ProtocolUnitResult]
    _unit_result_mapping: dict[ProtocolUnit, list[ProtocolUnitResult]]
    _result_unit_mapping: dict[ProtocolUnitResult, ProtocolUnit]

    def __init__(
        self,
        *,
        protocol_units: list[ProtocolUnit],
        protocol_unit_results: list[ProtocolUnitResult],
        transformation_key: GufeKey | None,
        extends_key: GufeKey | None = None,
        name: str | None = None,
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
        return {
            "name": self.name,
            "protocol_units": self._protocol_units,
            "protocol_unit_results": self._protocol_unit_results,
            "transformation_key": self._transformation_key,
            "extends_key": self._extends_key,
        }

    @classmethod
    def _from_dict(cls, dct: dict):
        return cls(**dct)

    @property
    def result_graph(self) -> nx.DiGraph:
        """
        DAG of `ProtocolUnitResult` nodes with edges denoting dependencies.

        Each edge is directed from a task towards its dependencies; for example,
        an edge between a production run and its equilibration would point
        towards the equilibration unit.
        """
        return self._result_graph

    @property
    def protocol_unit_results(self) -> list[ProtocolUnitResult]:
        """
        `ProtocolUnitResult`s for each `ProtocolUnit` used to compute this object.

        Results are given in DAG-dependency order. In this order, tasks are
        always listed after their dependencies.
        """
        return list(self._iterate_dag_order(self.result_graph))

    @property
    def protocol_unit_failures(self) -> list[ProtocolUnitFailure]:
        """A list of all failed units.

        Note
        ----
        These are returned in DAG order, with tasks listed after their
        dependencies.
        """
        # mypy can't figure out the types here, .ok() will ensure a certain type
        # https://mypy.readthedocs.io/en/stable/common_issues.html?highlight=cast#complex-type-tests
        return [r for r in self.protocol_unit_results if not r.ok()]  # type: ignore

    @property
    def protocol_unit_successes(self) -> list[ProtocolUnitResult]:
        """A list of only successful `ProtocolUnit` results.

        Note
        ----
        These are returned in DAG order, with tasks listed after their
        dependencies.
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
        MissingUnitResultError:
          if there are no results for that protocol unit
        ProtocolUnitFailureError:
          if there are only failures for that protocol unit
        """
        try:
            units = self._unit_result_mapping[protocol_unit]
        except KeyError:
            raise MissingUnitResultError(f"No such `protocol_unit`:{protocol_unit} present")
        else:
            for u in units:
                if u.ok():
                    return u
            else:
                raise ProtocolUnitFailureError(f"No success for `protocol_unit`:{protocol_unit} found")

    def unit_to_all_results(self, protocol_unit: ProtocolUnit) -> list[ProtocolUnitResult]:
        """Return all results (success and failure) for a given Unit.

        Returns
        -------
        results : list[ProtocolUnitResult]
          results for a given unit

        Raises
        ------
        MissingUnitResultError
          if no results present for a given unit
        """
        try:
            return self._unit_result_mapping[protocol_unit]
        except KeyError:
            raise MissingUnitResultError(f"No such `protocol_unit`:{protocol_unit} present")

    def result_to_unit(self, protocol_unit_result: ProtocolUnitResult) -> ProtocolUnit:
        try:
            return self._result_unit_mapping[protocol_unit_result]
        except KeyError:
            raise MissingUnitResultError(f"No such `protocol_unit_result`:{protocol_unit_result} present")

    def ok(self) -> bool:
        # ensure that for every protocol unit, there is an OK result object
        return all(any(pur.ok() for pur in self._unit_result_mapping[pu]) for pu in self._protocol_units)

    @property
    def terminal_protocol_unit_results(self) -> list[ProtocolUnitResult]:
        """Get ProtocolUnitResults that terminate the DAG.

        Returns
        -------
        list[ProtocolUnitResult]
          All ProtocolUnitResults which do not have a ProtocolUnitResult that
          follows on (depends) on them.
        """
        return [u for u in self._protocol_unit_results if not nx.ancestors(self._result_graph, u)]


class ProtocolDAG(GufeTokenizable, DAGMixin):
    """
    An executable directed acyclic graph (DAG) of :class:`ProtocolUnit` objects.

    A ``ProtocolDAG`` is composed of :class:`ProtocolUnit` objects as well as
    how they depend on each other. A single ``ProtocolDAG`` execution should
    yield sufficient information to calculate a free energy difference
    (though perhaps not converged) between two `ChemicalSystem` objects.

    A ``ProtocolDAG`` yields a ``ProtocolDAGResult`` when executed.

    Properties
    ----------
    name : str
        Optional identifier for this ``ProtocolDAGResult``.
    protocol_units : list[ProtocolUnit]
        ``ProtocolUnit`` s (given in DAG-dependency order) used to compute this
        ``ProtocolDAGResult``. Tasks are always listed after their dependencies.
    graph : nx.DiGraph
        Graph of ``ProtocolUnit`` s as nodes, with directed edges to each
        ``ProtocolUnit``'s dependencies.
    """

    def __init__(
        self,
        *,
        protocol_units: list[ProtocolUnit],
        transformation_key: GufeKey | None,
        extends_key: GufeKey | None = None,
        name: str | None = None,
    ):
        """Create a new `ProtocolDAG``

        Parameters
        ----------
        protocol_units : Iterable[ProtocolUnit]
            The `ProtocolUnit` s that make up this `ProtocolDAG`, with
            dependencies included as inputs.
        transformation_key : Optional[GufeKey]
            Key of the `Transformation` that this `ProtocolDAG` corresponds to, if
            applicable. This functions as a label for identifying the source of
            this `ProtocolDAG`. This label will be passed on to the
            `ProtocolDAGResult` resulting from execution of this `ProtocolDAG`.
        extends_key : Optional[GufeKey]
            Key of the `ProtocolDAGResult` that this `ProtocolDAG` extends from.
            This functions as a label for identifying the source of this
            `ProtocolDAG`. This label will be passed on to the
            `ProtocolDAGResult` resulting from execution of this `ProtocolDAG`.
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
        return {
            "name": self.name,
            "protocol_units": self.protocol_units,
            "transformation_key": self._transformation_key,
            "extends_key": self._extends_key,
        }

    @classmethod
    def _from_dict(cls, dct: dict):
        return cls(**dct)


def execute_DAG(
    protocoldag: ProtocolDAG,
    *,
    shared_basedir: Path,
    scratch_basedir: Path,
    stderr_basedir: Path,
    stdout_basedir: Path,
    keep_shared: bool = False,
    keep_scratch: bool = False,
    raise_error: bool = True,
    n_retries: int = 0,
) -> ProtocolDAGResult:
    """
    Locally execute a full :class:`ProtocolDAG` in serial and in-process.

    Parameters
    ----------
    protocoldag : ProtocolDAG
        The :class:``ProtocolDAG`` to execute.
    shared_basedir : Path
        Filesystem path to use for shared space that persists across whole DAG
        execution. Used by a `ProtocolUnit` to pass file contents to dependent
        class:``ProtocolUnit`` instances.
    scratch_basedir : Path
        Filesystem path to use for `ProtocolUnit` `scratch` space.
    stderr_basedir : Path
        Filesystem path to use for `ProtocolUnit` `stderr` archiving.
    stdout_basedir : Path
        Filesystem path to use for `ProtocolUnit` `stdout` archiving.
    keep_shared : bool
        If True, don't remove shared directories for `ProtocolUnit`s after
        the `ProtocolDAG` is executed.
    keep_scratch : bool
        If True, don't remove scratch directories for a `ProtocolUnit` after
        it is executed.
    raise_error : bool
        If True, raise an exception if a ProtocolUnit fails, default True
        if False, any exceptions will be stored as `ProtocolUnitFailure`
        objects inside the returned `ProtocolDAGResult`
    n_retries : int
        the number of times to attempt, default 0, i.e. try once and only once

    Returns
    -------
    ProtocolDAGResult
        The result of executing the `ProtocolDAG`.

    """
    if n_retries < 0:
        raise ValueError("Must give positive number of retries")

    # iterate in DAG order
    results: dict[GufeKey, ProtocolUnitResult] = {}
    all_results = []  # successes AND failures
    shared_paths = []
    for unit in protocoldag.protocol_units:
        # translate each `ProtocolUnit` in input into corresponding
        # `ProtocolUnitResult`
        inputs = _pu_to_pur(unit.inputs, results)

        attempt = 0
        while attempt <= n_retries:
            shared = shared_basedir / f"shared_{str(unit.key)}_attempt_{attempt}"
            shared_paths.append(shared)
            shared.mkdir()

            scratch = scratch_basedir / f"scratch_{str(unit.key)}_attempt_{attempt}"
            scratch.mkdir()

            stderr = stderr_basedir / f"stderr_{str(unit.key)}_attempt_{attempt}"
            stderr.mkdir()

            stdout = stdout_basedir / f"stdout_{str(unit.key)}_attempt_{attempt}"
            stdout.mkdir()

            context = Context(shared=shared, scratch=scratch, stderr=stderr, stdout=stdout)

            # execute
            result = unit.execute(context=context, raise_error=raise_error, **inputs)
            all_results.append(result)

            # clean up outputs
            shutil.rmtree(stderr)
            shutil.rmtree(stdout)

            if not keep_scratch:
                shutil.rmtree(scratch)

            if result.ok():
                # attach result to this `ProtocolUnit`
                results[unit.key] = result
                break
            attempt += 1

        if not result.ok():
            break

    if not keep_shared:
        for shared_path in shared_paths:
            shutil.rmtree(shared_path)

    return ProtocolDAGResult(
        name=protocoldag.name,
        protocol_units=protocoldag.protocol_units,
        protocol_unit_results=all_results,
        transformation_key=protocoldag.transformation_key,
        extends_key=protocoldag.extends_key,
    )


def _pu_to_pur(
    inputs: dict[str, Any] | list[Any] | ProtocolUnit,
    mapping: dict[GufeKey, ProtocolUnitResult],
):
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
