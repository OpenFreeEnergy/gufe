# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

"""Classes in this module must be subclassed by protocol authors.

"""

import abc
import warnings
from collections.abc import Iterable, Sized
from typing import Any, Optional, Union

from openff.units import Quantity

from ..chemicalsystem import ChemicalSystem
from ..mapping import ComponentMapping
from ..settings import Settings, SettingsBaseModel
from ..tokenization import GufeKey, GufeTokenizable
from .protocoldag import ProtocolDAG, ProtocolDAGResult
from .protocolunit import ProtocolUnit


class ProtocolResult(GufeTokenizable):
    """
    Container for all results for a single :class:`Transformation`.

    Contains a collection of :class:`ProtocolDAGResult` instances. This is an
    abstract base class; individual `Protocol` implementations should have a
    corresponding subclass of `ProtocolResult` implemented as well.

    The following methods should be implemented in any subclass:
    - `get_estimate`
    - `get_uncertainty`
    """

    def __init__(self, n_protocol_dag_results: int = 0, **data):
        self._data = data

        if not n_protocol_dag_results >= 0:
            raise ValueError("`n_protocol_dag_results` must be an integer greater than or equal to zero")

        self._n_protocol_dag_results = n_protocol_dag_results

    @classmethod
    def _defaults(cls):
        return {}

    def _to_dict(self):
        return {"n_protocol_dag_results": self.n_protocol_dag_results, "data": self.data}

    @classmethod
    def _from_dict(cls, dct: dict):
        # TODO: remove in gufe 2.0
        try:
            n_protocol_dag_results = dct["n_protocol_dag_results"]
        except KeyError:
            n_protocol_dag_results = 0
        return cls(n_protocol_dag_results=n_protocol_dag_results, **dct["data"])

    @property
    def n_protocol_dag_results(self) -> int:
        return self._n_protocol_dag_results

    @property
    def data(self) -> dict[str, Any]:
        """
        Aggregated data contents from multiple `ProtocolDAGResult` instances.

        The structure of this data is specific to the `Protocol` subclass each
        `ProtocolResult` subclass corresponds to.
        """
        return self._data

    @abc.abstractmethod
    def get_estimate(self) -> Quantity: ...

    @abc.abstractmethod
    def get_uncertainty(self) -> Quantity: ...


class Protocol(GufeTokenizable):
    """A method that via an alchemical transformation estimates free energy difference

    Takes a :class:`.Settings` object customised for this protocol on creation.
    This configures the protocol for repeated execution on (pairs of)
    :class:`ChemicalSystem` objects.

    This is an abstract base class; individual Protocol implementations
    should be subclasses of this class. The following methods should be
    implemented in any subclass:

    - `_create`
    - `_gather`
    - `_default_settings`

    Additionally, the `_settings_cls` must be set to the intended
    subclass of `SettingsBaseModel`. This attribute is validated
    during the instantiation of the Protocol.

    """

    _settings: Settings
    _settings_cls: type[SettingsBaseModel]
    result_cls: type[ProtocolResult]
    """Corresponding `ProtocolResult` subclass."""

    def __init__(self, settings: Settings):
        """
        Parameters
        ----------
        settings : Settings
          The parameters for this particular method.  This will be a specialised
          subclass for this particular Protocol.

        Note
        ----
        Once the Protocol object is created, the input Settings are frozen,
        so should be finalised before creating the Protocol instance.
        """

        if not hasattr(self.__class__, "_settings_cls"):
            raise NotImplementedError(
                f"class `{self.__class__.__qualname__}` must implement the `_settings_cls` attribute."
            )

        if not isinstance(settings, self._settings_cls):
            raise ValueError(
                f"`{self.__class__.__qualname__}` expected a `{self._settings_cls.__qualname__}` instance."
            )

        self._settings = settings.frozen_copy()

    @property
    def settings(self) -> Settings:
        """A read-only view of the settings for this ``Protocol`` instance."""
        return self._settings

    @classmethod
    def _defaults(cls):
        return {}

    def _to_dict(self):
        return {"settings": self.settings}

    @classmethod
    def _from_dict(cls, dct: dict):
        return cls(**dct)

    @classmethod
    @abc.abstractmethod
    def _default_settings(cls) -> Settings:
        """Method to override in custom `Protocol` subclasses.

        Gives a usable instance of ``Settings`` that function as
        reasonable defaults for this `Protocol` subclass.

        """
        raise NotImplementedError()

    @classmethod
    def default_settings(cls) -> Settings:
        """Get the default settings for this ``Protocol``.

        These represent the current best-practices for the use of this
        particular method. These can be modified and passed in as the only
        argument for creating a new ``Protocol`` instance.
        """
        return cls._default_settings()

    @abc.abstractmethod
    def _create(
        self,
        stateA: ChemicalSystem,
        stateB: ChemicalSystem,
        mapping: ComponentMapping | list[ComponentMapping] | None,
        extends: ProtocolDAGResult | None = None,
    ) -> list[ProtocolUnit]:
        """Method to override in custom :class:`Protocol` subclasses.

        This method should take two `ChemicalSystem`s, and optionally one or
         more ``ComponentMapping`` objects, and prepare a collection of
        ``ProtocolUnit`` instances that when executed in order give sufficient
        information to estimate the free energy difference between those two
        `ChemicalSystem`s.

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
        *,
        stateA: ChemicalSystem,
        stateB: ChemicalSystem,
        mapping: ComponentMapping | list[ComponentMapping] | dict[str, ComponentMapping] | None,
        extends: ProtocolDAGResult | None = None,
        name: str | None = None,
        transformation_key: GufeKey | None = None,
    ) -> ProtocolDAG:
        r"""Prepare a `ProtocolDAG` with all information required for execution.

        A :class:`.ProtocolDAG` is composed of :class:`.ProtocolUnit` \s, with
        dependencies established between them. These form a directed, acyclic
        graph, and each `ProtocolUnit` can be executed once its dependencies have
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
        mapping : Optional[Union[ComponentMapping, list[ComponentMapping]]]
            Mappings of e.g. atoms between a labelled component in the
            stateA and stateB `ChemicalSystem` .
        extends : Optional[ProtocolDAGResult]
            If provided, then the `ProtocolDAG` produced will start from the
            end state of the given `ProtocolDAGResult`. This allows for
            extension from a previously-run `ProtocolDAG`.
        name : Optional[str]
            A user supplied identifier for the resulting DAG
        transformation_key : Optional[GufeKey]
            Key of the `Transformation` that this `Protocol` corresponds to, if
            applicable. This will be used to label the resulting `ProtocolDAG`,
            and can be used for identifying its source. This label will be
            passed on to the `ProtocolDAGResult` resulting from execution of
            this `ProtocolDAG`.

        Returns
        -------
        ProtocolDAG
            A directed, acyclic graph that can be executed by a `Scheduler`.
        """
        if isinstance(mapping, dict):
            warnings.warn(
                ("mapping input as a dict is deprecated, " "instead use either a single Mapping or list"),
                DeprecationWarning,
            )
            mapping = list(mapping.values())

        return ProtocolDAG(
            name=name,
            protocol_units=self._create(
                stateA=stateA,
                stateB=stateB,
                mapping=mapping,
                extends=extends,
            ),
            transformation_key=transformation_key,
            extends_key=extends.key if extends is not None else None,
        )

    def gather(self, protocol_dag_results: Iterable[ProtocolDAGResult]) -> ProtocolResult:
        """Gather multiple ProtocolDAGResults into a single ProtocolResult.

        Parameters
        ----------
        protocol_dag_results : Iterable[ProtocolDAGResult]
            The `ProtocolDAGResult` objects to assemble aggregate quantities
            from.

        Returns
        -------
        ProtocolResult
            Aggregated results from many `ProtocolDAGResult`s from a given `Protocol`.
        """
        # Iterable does not implement __len__ and makes no guarantees that
        # protocol_dag_results is finite, checking both in method signature
        # doesn't appear possible, explicitly check for __len__ through the
        # Sized type
        if not isinstance(protocol_dag_results, Sized):
            raise ValueError("`protocol_dag_results` must implement `__len__`")
        return self.result_cls(n_protocol_dag_results=len(protocol_dag_results), **self._gather(protocol_dag_results))

    @abc.abstractmethod
    def _gather(self, protocol_dag_results: Iterable[ProtocolDAGResult]) -> dict[str, Any]:
        """Method to override in custom Protocol subclasses.

        This method should take any number of ``ProtocolDAGResult``s produced
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
