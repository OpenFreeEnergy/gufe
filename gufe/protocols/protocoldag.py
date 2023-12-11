# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from copy import copy
from collections import defaultdict
import os
from typing import Iterable, Optional, Union, Any
from os import PathLike
from pathlib import Path
import shutil

import networkx as nx

from ..tokenization import GufeTokenizable, GufeKey
from .protocolunit import (
    ProtocolUnit, ProtocolUnitResult, ProtocolUnitFailure, Context
)

from ..storage.storagemanager import StorageManager
from ..storage.externalresource.filestorage import FileStorage
from ..storage.externalresource.base import ExternalStorage

import logging
_logger = logging.getLogger(__name__)


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
    def name(self) -> Optional[str]:
        """Optional identifier."""
        return self._name

    @property
    def graph(self) -> nx.DiGraph:
        """DAG of `ProtocolUnit` nodes with edges denoting dependencies."""
        return self._graph

    @property
    def protocol_units(self) -> list[ProtocolUnit]:
        """
        List of `ProtocolUnit`s given in DAG-dependency order.

        DAG-dependency order guarantees that any task is listed after all of its
        dependencies.
        """
        return list(self._iterate_dag_order(self._graph))

    @property
    def transformation_key(self) -> Union[GufeKey, None]:
        """
        The `GufeKey` of the `Transformation` this object performs.

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
    """
    An executable directed acyclic graph (DAG) of :class:`ProtocolUnit` objects.

    A ``ProtocolDAG`` is composed of :class:`ProtocolUnit` objects as well as
    how they depend on each other. A single `ProtocolDAG` execution should
    yield sufficient information to calculate a free energy difference
    (though perhaps not converged) between two `ChemicalSystem` objects.
    
    A `ProtocolDAG` yields a `ProtocolDAGResult` when executed.

    Attributes
    ----------
    name : str
        Optional identifier for this `ProtocolDAGResult`.
    protocol_units : list[ProtocolUnit]
        `ProtocolUnit`s (given in DAG-dependency order) used to compute this
        `ProtocolDAGResult`. Tasks are always listed after their dependencies.
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


class ReproduceOldBehaviorStorageManager(StorageManager):
    # Default behavior has scratch at {dag_label}/scratch/{unit_label} and
    # shared at {dag_label}/{unit_label}. This little class makes changes
    # that get us back to the original behavior of this class: scratch at
    # {dag_label}/scratch_{unit_label} and shared at
    # {dag_label}/shared_{unit_label}.
    def _scratch_loc(self, dag_label, unit_label, attempt):
        return (
            self.scratch_root
            / f"{dag_label}/scratch_{unit_label}_attempt_{attempt}"
        )

    def make_label(self, dag_label, unit_label, attempt):
        return f"{dag_label}/shared_{unit_label}_attempt_{attempt}"

    @classmethod
    def from_old_args(
        cls,
        shared_basedir: PathLike,
        scratch_basedir: PathLike, *,
        keep_shared: bool = False,
        keep_scratch: bool = False,
    ):
        """
        Create an new storage manager based on the old execute_DAG args.

        Parameters
        ----------
        shared_basedir : Path
            Filesystem path to use for shared space that persists across whole DAG
            execution. Used by a `ProtocolUnit` to pass file contents to dependent
            class:``ProtocolUnit`` instances.
        scratch_basedir : Path
            Filesystem path to use for `ProtocolUnit` `scratch` space.
        keep_shared : bool
            If True, don't remove shared directories for `ProtocolUnit`s after
            the `ProtocolDAG` is executed.
        keep_scratch : bool
            If True, don't remove scratch directories for a `ProtocolUnit` after
            it is executed.
        """
        # doing this here makes it easier to test than putting in
        # execute_DAG
        shared_basedir = Path(shared_basedir)
        shared = FileStorage(shared_basedir.parent)
        storage_manager = cls(
            scratch_root=scratch_basedir,
            shared_root=shared,
            permanent_root=shared,
            keep_scratch=keep_scratch,
            keep_shared=keep_shared,
            keep_staging=True,
            keep_empty_dirs=True,
            staging=Path(""),  # use the actual directories as the staging
        )
        return storage_manager


def execute_DAG(protocoldag: ProtocolDAG, *,
                shared_basedir: PathLike,
                scratch_basedir: PathLike,
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
    # the directory given as shared_root is actually the directory for this
    # DAG; the "shared_root" for the storage manager is the parent. We'll
    # force permanent to be the same.
    storage_manager = ReproduceOldBehaviorStorageManager.from_old_args(
        shared_basedir=shared_basedir,
        scratch_basedir=scratch_basedir,
        keep_shared=keep_shared,
        keep_scratch=keep_scratch
    )
    dag_label = shared_basedir.name
    return new_execute_DAG(protocoldag, dag_label, storage_manager,
                           raise_error, n_retries)


def new_execute_DAG(  # TODO: this is a terrible name
    protocoldag: ProtocolDAG,
    dag_label: str,
    storage_manager: StorageManager,
    raise_error: bool = False,
    n_retries: int = 0
) -> ProtocolDAGResult:
    """
    Locally execute a full :class:`ProtocolDAG` in serial and in-process.

    Alternate input signature to generalize execute_DAG

    Parameters
    ----------
    protocoldag : ProtocolDAG
        The :class:``ProtocolDAG`` to execute.
    dag_label : str
        Label to use for the DAG
    storage_manager : StorageManager
        The :class:`.StorageManager` to handle storing files.
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
    # this simplifies setup of execute_DAG by allowing you to directly
    # provide the storage_manager; the extra options in the old one just
    # configure the storage_manager
    if n_retries < 0:
        raise ValueError("Must give positive number of retries")

    # iterate in DAG order
    results: dict[GufeKey, ProtocolUnitResult] = {}
    all_results = []  # successes AND failures

    with storage_manager.running_dag(dag_label) as dag_ctx:
        for unit in protocoldag.protocol_units:
            # import pdb; pdb.set_trace()
            attempt = 0
            while attempt <= n_retries:
                # translate each `ProtocolUnit` in input into corresponding
                # `ProtocolUnitResult`
                inputs = _pu_to_pur(unit.inputs, results)

                label = storage_manager.make_label(dag_label, unit.key,
                                                   attempt=attempt)
                with dag_ctx.running_unit(
                    dag_label, unit.key, attempt=attempt
                ) as context:
                    _logger.info("Starting unit {label}")
                    _logger.info(context)
                    result = unit.execute(
                            context=context,
                            raise_error=raise_error,
                            **inputs)
                    all_results.append(result)

                if result.ok():
                    # attach result to this `ProtocolUnit`
                    results[unit.key] = result
                    break
                attempt += 1

            if not result.ok():
                break

    return ProtocolDAGResult(
            name=protocoldag.name, 
            protocol_units=protocoldag.protocol_units, 
            protocol_unit_results=all_results,
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

