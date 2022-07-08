# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

"""Classes in this module must be subclassed by protocol authors.

"""

import abc
from typing import Optional, Iterable, Any, Dict

from openff.toolkit.utils.serialization import Serializable
import networkx as nx

from ..chemicalsystem import ChemicalSystem
from ..mapping import Mapping

from .protocoldag import ProtocolDAG
from .results import ProtocolDAGResult


class ProtocolResult(Serializable, abc.ABC):
    """Container for all `ProtocolDAGResult`s for a given `Transformation`."""

    def __init__(self, data, **kwargs):
        self._data = data

    @property
    def data(self):
        return self._data

    def to_dict(self) -> dict:
        ...

    @classmethod
    def from_dict(cls, d: dict):
        ...

    @abc.abstractmethod
    def get_estimate(self):
        ...

    @abc.abstractmethod
    def get_uncertainty(self):
        ...

    @abc.abstractmethod
    def get_rate_of_convergence(self):
        ...


class Protocol(Serializable, abc.ABC):
    """A protocol that implements an alchemical transformation.

    Takes a `ProtocolSettings` object specific to the protocol on init.
    This configures the protocol for repeated execution on `ChemicalSystem`s.

    Attributes
    ----------
    settings : ProtocolSettings

    """

    _results_cls = ProtocolResult

    def __init__(self, settings: "ProtocolSettings" = None):  # type: ignore
        """ """
        self._settings = settings

    @property
    def settings(self):
        return self._settings

    def to_dict(self) -> dict:
        ...

    @classmethod
    def from_dict(cls, d: dict):
        ...

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        if self._settings != other.settings:
            return False

        return True

    def __hash__(self):
        return hash((self.__class__.__name__, self._settings))

    @classmethod
    @abc.abstractmethod
    def get_default_settings(cls):
        """Get the default settings for this protocol.

        These can be modified and passed back in to the class init.

        """
        ...

    @abc.abstractmethod
    def _create(
        self,
        stateA: ChemicalSystem,
        stateB: ChemicalSystem,
        mapping: Optional[Mapping] = None,
        extend_from: Optional[ProtocolDAGResult] = None,
    ) -> nx.DiGraph:
        ...

    def create(
        self,
        stateA: ChemicalSystem,
        stateB: ChemicalSystem,
        mapping: Optional[Mapping] = None,
        extend_from: Optional[ProtocolDAGResult] = None,
        name: str = None,
    ) -> ProtocolDAG:
        """Prepare a `ProtocolDAG` with all information required for execution.

        A `ProtocolDAG` is composed of `ProtocolUnit`s, with dependencies
        established between them. These form a directed, acyclic graph,
        and each `ProtocolUnit` can be executed once its dependencies have
        completed.

        A `ProtocolDAG` can be passed to a `Scheduler` for execution on its
        resources. A `ProtocolResult` can be retrieved from the `Scheduler`
        upon completion of all `ProtocolUnit`s in the `ProtocolDAG`.

        Parameters
        ----------
        stateA : ChemicalSystem
            The starting `ChemicalSystem` for the transformation.
        stateB : ChemicalSystem
            The ending `ChemicalSystem` for the transformation.
        mapping : Optional[Mapping]
            Mapping of e.g. atoms between the `stateA` and `stateB`
            `ChemicalSystem`s.
        extend_from : Optional[ProtocolDAGResult]
            If provided, then the `ProtocolDAG` produced will start from the
            end state of the given `ProtocolDAGResult`. This allows for
            extension from a previously-run `ProtocolDAG`.
        name : Optional[str]
            A user supplied identifier for the resulting DAG

        Returns
        -------
        ProtocolDAG
            A directed, acyclic graph that can be executed by a `Scheduler`.

        """
        return ProtocolDAG(
            name=name,
            protocol_units=self._create(
                stateA=stateA,
                stateB=stateB,
                mapping=mapping,
                extend_from=extend_from,
            ),
        )

    def gather(
        self, protocol_dag_results: Iterable[ProtocolDAGResult]
    ) -> ProtocolResult:
        """Gather multiple `ProtocolDAGResult`s into a single `ProtocolResult`.

        Parameters
        ----------
        protocol_dag_results : Iterable[ProtocolDAGResult]
            The `ProtocolDAGResult`s to assemble aggregate quantities from.

        Returns
        -------
        ProtocolResult
            Aggregated results from many `ProtocolDAGResult`s from a given `Protocol`.

        """
        return self.result_cls(**self._gather(protocol_dag_results))

    @abc.abstractmethod
    def _gather(
        self, protocol_dag_results: Iterable[ProtocolDAGResult]
    ) -> Dict[str, Any]:
        ...
