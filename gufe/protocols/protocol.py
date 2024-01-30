# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

"""Classes in this module must be subclassed by protocol authors.

"""

import abc
from typing import Optional, Iterable, Any, Union
from openff.units import Quantity
import functools
import json

from ..settings import Settings
from ..tokenization import GufeTokenizable, GufeKey, JSON_HANDLER
from ..custom_json import JSONSerializerDeserializer, JSONCodec
from ..custom_codecs import default_from_dict, inherited_is_my_dict
from gufe.settings.models import SettingsBaseModel
from ..chemicalsystem import ChemicalSystem
from ..mapping import ComponentMapping

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

    def __init__(self, **data):
        self._data = data

    @classmethod
    def _defaults(cls):
        return {}

    def _to_dict(self):
        return {'data': self.data}

    @classmethod
    def _from_dict(cls, dct: dict):
        return cls(**dct['data'])

    @property
    def data(self) -> dict[str, Any]:
        """
        Aggregated data contents from multiple `ProtocolDAGResult` instances.

        The structure of this data is specific to the `Protocol` subclass each
        `ProtocolResult` subclass corresponds to.
        """
        return self._data

    @abc.abstractmethod
    def get_estimate(self) -> Quantity:
        ...

    @abc.abstractmethod
    def get_uncertainty(self) -> Quantity:
        ...


class Protocol(GufeTokenizable):
    """A protocol that implements an alchemical transformation.

    Takes a ``Settings`` object customised for this protocol on init.
    This configures the protocol for repeated execution on (pairs of)
    :class:`ChemicalSystem` objects.

    This is an abstract base class; individual Protocol implementations
    should be subclasses of this class. The following methods should be
    implemented in any subclass:

    - `_create`
    - `_gather`
    - `_default_settings`
    """
    _settings: Settings
    result_cls: type[ProtocolResult]
    """Corresponding `ProtocolResult` subclass."""

    def __init__(self, settings: Settings):
        """Create a new ``Protocol`` instance.

        Parameters
        ----------
        settings : Settings
            The full settings for this ``Protocol`` instance.
        """
        self._settings = settings

    @property
    def settings(self) -> Settings:
        """The full settings for this ``Protocol`` instance."""
        return self._settings

    @classmethod
    def _defaults(cls):
        return {}

    def _to_dict(self):
        return {'settings': self.settings}

    @classmethod
    def _from_dict(cls, dct: dict):
        # This is an ugly band-aid on the problem, but it enables settings
        # migration
        settings_codec = JSONCodec(
            cls=SettingsBaseModel,
            to_dict=lambda obj: {field: getattr(obj, field)
                                 for field in obj.__fields__},
            from_dict=default_from_dict,
            is_my_dict=functools.partial(inherited_is_my_dict,
                                         cls=SettingsBaseModel),
        )
        # by putting ours first, it gets hit before the one that shows up
        # later
        handler = JSONSerializerDeserializer([settings_codec] +
                                             JSON_HANDLER.codecs)
        json_str = json.dumps(dct, cls=JSON_HANDLER.encoder)
        updated_dct = json.loads(json_str, cls=handler.decoder)
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
        extends: Optional[ProtocolDAGResult] = None,
    ) -> list[ProtocolUnit]:
        """Method to override in custom :class:`Protocol` subclasses.

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
        *,
        stateA: ChemicalSystem,
        stateB: ChemicalSystem,
        mapping: Union[dict[str, ComponentMapping], None],
        extends: Optional[ProtocolDAGResult] = None,
        name: Optional[str] = None,
        transformation_key: Optional[GufeKey] = None
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
        return ProtocolDAG(
            name=name,
            protocol_units=self._create(
                stateA=stateA,
                stateB=stateB,
                mapping=mapping,
                extends=extends,
            ),
            transformation_key=transformation_key,
            extends_key=extends.key if extends is not None else None
        )

    def gather(
        self, protocol_dag_results: Iterable[ProtocolDAGResult]
    ) -> ProtocolResult:
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
        return self.result_cls(**self._gather(protocol_dag_results))

    @abc.abstractmethod
    def _gather(
        self, protocol_dag_results: Iterable[ProtocolDAGResult]
    ) -> dict[str, Any]:
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
