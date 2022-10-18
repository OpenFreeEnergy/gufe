# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

"""Classes in this module must be subclassed by protocol authors.

"""

import abc
from typing import Optional, Iterable, Any, Dict, List

import networkx as nx

from ..tokenization import GufeTokenizable
from ..chemicalsystem import ChemicalSystem
from ..mapping import ComponentMapping

from .protocoldag import ProtocolDAG, ProtocolDAGResult
from .protocolunit import ProtocolUnit


class ProtocolResult(GufeTokenizable):
    """Container for all `ProtocolDAGResult`s for a given `Transformation`.

    This is an abstract base class; individual `Protocol` implementations
    should have a corresponding subclass of `ProtocolResult` implemented as
    well. 

    The following methods should be implemented in any subclass:
    - `get_estimate`
    - `get_uncertainty`
    - `get_rate_of_convergence`

    Attributes
    ----------
    data : Dict[str,Any]
        Aggregated data contents from multiple `ProtocolDAGResult`s.
        The structure of this data is specific to the `Protocol` subclass each
        `ProtocolResult` subclass corresponds to.

    """

    def __init__(self, **data):
        self._data = data

    def _defaults(self):
        return {}

    def _to_dict(self):
        return {'data': self.data}

    @classmethod
    def _from_dict(cls, dct: Dict):
        return cls(**dct)

    @property
    def data(self):
        return self._data

    @abc.abstractmethod
    def get_estimate(self):
        ...

    @abc.abstractmethod
    def get_uncertainty(self):
        ...

    @abc.abstractmethod
    def get_rate_of_convergence(self):
        ...


class Protocol(GufeTokenizable):
    """A protocol that implements an alchemical transformation.

    Takes a `ProtocolSettings` object specific to the protocol on init.
    This configures the protocol for repeated execution on `ChemicalSystem`s.

    This is an abstract base class; individual `Protocol` implementations
    should be subclasses of this class. The following methods should be
    implemented in any subclass:
    - `_create`
    - `_gather`
    - `_default_settings`

    Attributes
    ----------
    settings : ProtocolSettings
        The full settings for this `Protocol` instance.
    result_cls : type[ProtocolResult]
        Correponding `ProtocolResult` subclass t

    """

    result_cls: type[ProtocolResult]

    def __init__(self, settings: "ProtocolSettings" = None):  # type: ignore
        """Create a new `Protocol` instance.

        Parameters
        ----------
        settings : ProtocolSettings
            The full settings for this `Protocol` instance.

        """
        self._settings = settings

    @property
    def settings(self):
        return self._settings

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        if self._settings != other.settings:
            return False

        return True

    def __hash__(self):
        return hash((self.__class__.__name__, self._settings))

    def _defaults(self):
        return {}

    def _to_dict(self):
        return {'settings': self.settings}

    @classmethod
    def _from_dict(cls, dct: Dict):
        return cls(**dct)

    @classmethod
    @abc.abstractmethod
    def _default_settings(cls) -> "ProtocolSettings":     # type: ignore
        """Method to override in custom `Protocol` subclasses.

        Gives a usable instance of `ProtocolSettings` that function as
        reasonable defaults for this `Protocol` subclass.

        """
        raise NotImplementedError()

    @classmethod
    def default_settings(cls) -> "ProtocolSettings":      # type: ignore
        """Get the default settings for this `Protocol`.

        These can be modified and passed in as the `settings` for a new
        `Protocol` instance.

        """
        return cls._default_settings()

    @abc.abstractmethod
    def _create(
        self,
        stateA: ChemicalSystem,
        stateB: ChemicalSystem,
        mapping: Optional[dict[str, ComponentMapping]] = None,
        extend_from: Optional[ProtocolDAGResult] = None,
    ) -> List[ProtocolUnit]:
        """Method to override in custom `Protocol` subclasses.

        This method should take two `ChemicalSystem`s, and optionally a
        dict mapping string to ``ComponentMapping``, and prepare a collection of ``ProtocolUnit`` instances
        that when executed in order give sufficient information to estimate the
        free energy difference between those two `ChemicalSystem`s.

        This method should return a list of `ProtocolUnit` instances.
        For an instance in which another `ProtocolUnit` is given as a parameter
        on init, the former instance will be executed only *after* the latter,
        with the latter's `ProtocolUnitResult` provided as corresponding input
        to the former's `execute` method.

        In this way, the list of `ProtocolUnit`s returned by this method give an
        implicit dependency DAG (directed, acyclic graph).

        This method can optionally support extension from an existing
        `ProtocolDAGResult` by extracting the information needed to use it as a
        starting point for any simulations required by the `Protocol`.

        See also
        --------
        :meth:`Protocol.create`

        """
        ...

    def create(
        self,
        stateA: ChemicalSystem,
        stateB: ChemicalSystem,
        mapping: Optional[dict[str, ComponentMapping]] = None,
        extend_from: Optional[ProtocolDAGResult] = None,
        name: str = None,
    ) -> ProtocolDAG:
        """Prepare a `ProtocolDAG` with all information required for execution.

        A `ProtocolDAG` is composed of `ProtocolUnit`s, with dependencies
        established between them. These form a directed, acyclic graph,
        and each `ProtocolUnit` can be executed once its dependencies have
        completed.

        A `ProtocolDAG` can be passed to a `Scheduler` for execution on its
        resources. A `ProtocolDAGResult` can be retrieved from the `Scheduler`
        upon completion of all `ProtocolUnit`s in the `ProtocolDAG`.

        Parameters
        ----------
        stateA : ChemicalSystem
            The starting `ChemicalSystem` for the transformation.
        stateB : ChemicalSystem
            The ending `ChemicalSystem` for the transformation.
        mapping : Optional[dict[str, ComponentMapping]]
            Mappings of e.g. atoms between a labelled component in the
             `stateA` and `stateB` `ChemicalSystem`s.
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
        """Method to override in custom `Protocol` subclasses.

        This method should take any number of `ProtocolDAGResult`s produced
        by this `Protocol` for a pair of `ChemicalSystem`s and extract the
        quantities necessary to calculate aggregated estimates of the free
        energy difference between them.

        This method returns a dict, which becomes the `data` attribute of this
        `Protocol` subclass's corresponding `ProtocolResult` subclass.

        See also
        --------
        :meth:`Protocol.gather`

        """
        ...
