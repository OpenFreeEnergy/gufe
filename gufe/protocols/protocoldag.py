# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any
from json import JSONDecodeError

import networkx as nx

from ..tokenization import GufeKey, GufeTokenizable
from .errors import MissingUnitResultError, ProtocolDAGError, ProtocolUnitFailureError
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
        """List of ``ProtocolUnit`` s given in DAG-dependency order.

        DAG-dependency order guarantees that any task is listed after all of its
        dependencies.
        """
        return list(self._iterate_dag_order(self._graph))

    @property
    def transformation_key(self) -> GufeKey | None:
        """The ``GufeKey`` of the ``Transformation`` this object performs.

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
    """Result for a single execution of an entire :class:`ProtocolDAG`.

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

        # build mappings from protocol units to results
        unit_key_to_unit = {unit.key: unit for unit in protocol_units}
        unit_key_to_results = defaultdict(list)
        self._result_unit_mapping = dict()
        self._unit_result_mapping = dict()

        for result in protocol_unit_results:
            unit_key_to_results[result.source_key].append(result)
            self._result_unit_mapping[result] = unit_key_to_unit[result.source_key]

        for unit in protocol_units:
            self._unit_result_mapping[unit] = unit_key_to_results.get(unit.key, list())

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
        """DAG of ``ProtocolUnitResult`` nodes with edges denoting dependencies.

        Each edge is directed from a task towards its dependencies; for example,
        an edge between a production run and its equilibration would point
        towards the equilibration unit.
        """
        return self._result_graph

    @property
    def protocol_unit_results(self) -> list[ProtocolUnitResult]:
        r"""`ProtocolUnitResult`\s for each `ProtocolUnit` used to compute this object.

        Results are given in DAG-dependency order. In this order, tasks are
        always listed after their dependencies.
        """
        return list(self._iterate_dag_order(self.result_graph))

    @property
    def protocol_unit_failures(self) -> list[ProtocolUnitFailure]:
        r"""A list of all failed ``ProtocolUnit``\s.

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
        """Return the successful ``ProtocolUnitResult`` for a given ``ProtocolUnit``,
        if it exists.

        Returns
        -------
        success : ProtocolUnitResult
          the successful ``ProtocolUnitResult`` for this ``ProtocolUnit``

        Raises
        ------
        MissingUnitResultError:
          if there are no results for that ``ProtocolUnit``
        ProtocolUnitFailureError:
          if there are only failures for that ``ProtocolUnit``
        """
        try:
            unit_results = self._unit_result_mapping[protocol_unit]
        except KeyError:
            raise MissingUnitResultError(f"No such `protocol_unit`:{protocol_unit} present")
        else:
            for u in unit_results:
                if u.ok():
                    return u
            else:
                raise ProtocolUnitFailureError(f"No success for `protocol_unit`:{protocol_unit} found")

    def unit_to_all_results(self, protocol_unit: ProtocolUnit) -> list[ProtocolUnitResult]:
        r"""Return all ``ProtocolUnitResult``\s (success and failure) for a given ``ProtocolUnit``.

        Returns
        -------
        results : list[ProtocolUnitResult]
          ``ProtocolUnitResult``\s for the given ``ProtocolUnit``

        Raises
        ------
        MissingUnitResultError
          if no ``ProtocolUnitResult``\s present for the given ``ProtocolUnit``
        """
        try:
            return self._unit_result_mapping[protocol_unit]
        except KeyError:
            raise MissingUnitResultError(f"No such `protocol_unit`:{protocol_unit} present")

    def result_to_unit(self, protocol_unit_result: ProtocolUnitResult) -> ProtocolUnit:
        """Return the ``ProtocolUnit`` corresponding to the given ``ProtocolUnitResult``.

        Returns
        -------
        ProtocolUnit
            the ``ProtocolUnit`` corresponding to the given ``ProtocolUnitResult``

        Raises
        ------
        MissingUnitResultError
          if the given ``ProtocolUnitResult`` isn't present
        """
        try:
            return self._result_unit_mapping[protocol_unit_result]
        except KeyError:
            raise MissingUnitResultError(f"No such `protocol_unit_result`:{protocol_unit_result} present")

    def ok(self) -> bool:
        """Check if this ``ProtocolDAGResult`` succeeded or failed.

        Returns ``True`` if there is at least one successful ``ProtocolUnitResult`` for each ``ProtocolUnit``,
        and ``False`` otherwise.

        """
        # ensure that for every protocol unit, there is an OK result object
        return all(any(pur.ok() for pur in self._unit_result_mapping[pu]) for pu in self._protocol_units)

    @property
    def terminal_protocol_unit_results(self) -> list[ProtocolUnitResult]:
        r"""Get ``ProtocolUnitResult``\s that terminate the DAG.

        Returns
        -------
        list[ProtocolUnitResult]
          All ``ProtocolUnitResult``\s which do not have a ``ProtocolUnitResult``
          that follows on (depends) on them.
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

        # validate that all units in the graph were explicitly provided
        provided_units = set(protocol_units)
        graph_units = set(self._graph.nodes)
        missing_units = graph_units - provided_units
        if missing_units:
            raise ProtocolDAGError(
                f"ProtocolDAG contains units that were not explicitly provided: "
                f"{missing_units}. All units must be passed explicitly via `protocol_units`, "
                f"even if they are dependencies of other units."
            )

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
    unitresults_basedir: Path | None = None,
    stderr_basedir: Path | None = None,
    stdout_basedir: Path | None = None,
    keep_shared: bool = False,
    keep_scratch: bool = False,
    keep_unitresults: bool = False,
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
    unitresults_basedir : Path | None = None
        Filesystem path to use for `ProtocolUnitResult` archiving during
        execution. If ``None`` (default), results will not be archived
        and it will not be able to resume DAG execution from the last
        succesfully finished `ProtocolUnit`.
    stderr_basedir : Path | None
        Filesystem path to use for `ProtocolUnit` `stderr` archiving.
    stdout_basedir : Path | None
        Filesystem path to use for `ProtocolUnit` `stdout` archiving.
    keep_shared : bool
        If True, don't remove shared directories for `ProtocolUnit`s after
        the `ProtocolDAG` is executed.
    keep_scratch : bool
        If True, don't remove scratch directories for a `ProtocolUnit` after
        it is executed.
    keep_unitresults : bool
        If True, don't remove the unitresults directory which contains
        the serialized `ProtocolUnitResult` for all executed `ProtocolUnit`.
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

    # handle results & optionally archiving
    results: dict[GufeKey, ProtocolUnitResult] = {}

    if unitresults_basedir is not None:
        unitresults_path = unitresults_basedir / f"unitresults_{str(protocoldag.key)}"
        unitresults_path.mkdir(exist_ok=True, parents=True)

        for file in unitresults_path.rglob("*.json"):
            try:
                unit_result = ProtocolUnitResult.from_json(file)
            except JSONDecodeError:
                pass
            else:
                # Is source key stable enough?
                # We probably don't want to resume if gufe stability has changed
                results[unit_result.source_key] = unit_result
        

    # iterate in DAG order
    all_results = []  # successes AND failures
    shared_paths = []
    for unit in protocoldag.protocol_units:
        # If we already have results, skip execution
        if unit.key in results:
            continue

        # translate each `ProtocolUnit` in input into corresponding
        # `ProtocolUnitResult`
        inputs = _pu_to_pur(unit.inputs, results)

        attempt = 0
        while attempt <= n_retries:
            shared = shared_basedir / f"shared_{str(unit.key)}_attempt_{attempt}"
            shared_paths.append(shared)
            shared.mkdir(exist_ok=True)

            scratch = scratch_basedir / f"scratch_{str(unit.key)}_attempt_{attempt}"
            scratch.mkdir(exist_ok=True)

            stderr = None
            if stderr_basedir:
                stderr = stderr_basedir / f"stderr_{str(unit.key)}_attempt_{attempt}"
                stderr.mkdir(exist_ok=True)

            stdout = None
            if stdout_basedir:
                stdout = stdout_basedir / f"stdout_{str(unit.key)}_attempt_{attempt}"
                stdout.mkdir(exist_ok=True)

            context = Context(shared=shared, scratch=scratch, stderr=stderr, stdout=stdout)

            # execute
            result = unit.execute(context=context, raise_error=raise_error, **inputs)
            all_results.append(result)

            # clean up outputs
            if stderr:
                shutil.rmtree(stderr)
            if stdout:
                shutil.rmtree(stdout)

            if not keep_scratch:
                shutil.rmtree(scratch)

            if result.ok():
                # attach result to this `ProtocolUnit`
                results[unit.key] = result

                # Serialize results if requested
                if unitresults_basedir is not None:
                    result.to_json(unitresults_path / f"{str(result.key)}.json")

                break
            attempt += 1

        if not result.ok():
            break

    if not keep_shared:
        for shared_path in shared_paths:
            shutil.rmtree(shared_path)

    if not keep_unitresults and unitresults_basedir is not None:
        shutil.rmtree(unitresults_path)

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
