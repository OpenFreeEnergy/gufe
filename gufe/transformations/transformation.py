# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
import json
import warnings
from collections.abc import Iterable
from typing import Optional, Union

from ..chemicalsystem import ChemicalSystem
from ..mapping import ComponentMapping
from ..protocols import Protocol, ProtocolDAG, ProtocolDAGResult, ProtocolResult
from ..tokenization import JSON_HANDLER, GufeTokenizable
from ..utils import ensure_filelike


class TransformationBase(GufeTokenizable):
    def __init__(
        self,
        protocol: Protocol,
        name: Optional[str] = None,
    ):
        """Transformation base class.

        Parameters
        ----------
        protocol : Protocol
            The sampling method to use for the transformation.
        name : str, optional
            A human-readable name for this transformation.
    
        """
        self._protocol = protocol
        self._name = name

    @classmethod
    def _defaults(cls):
        return super()._defaults()

    @property
    def name(self) -> Optional[str]:
        """Optional identifier for the transformation; used as part of its hash.

        Set this to a unique value if adding multiple, otherwise identical
        transformations to the same :class:`.AlchemicalNetwork` to avoid
        deduplication.
        """
        return self._name

    @classmethod
    def _from_dict(cls, d: dict):
        return cls(**d)

    @property
    @abc.abstractmethod
    def stateA(self) -> ChemicalSystem:
        """The starting :class:`.ChemicalSystem` for the transformation."""
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def stateB(self) -> ChemicalSystem:
        """The ending :class:`.ChemicalSystem` for the transformation."""
        raise NotImplementedError

    @abc.abstractmethod
    def create(
        self,
        *,
        extends: Optional[ProtocolDAGResult] = None,
        name: Optional[str] = None,
    ) -> ProtocolDAG:
        """
        Returns a :class:`.ProtocolDAG` executing this ``Transformation.protocol``.
        """
        raise NotImplementedError

    def gather(self, protocol_dag_results: Iterable[ProtocolDAGResult]) -> ProtocolResult:
        """Gather multiple :class:`.ProtocolDAGResult` \s into a single
        :class:`.ProtocolResult`.

        Parameters
        ----------
        protocol_dag_results : Iterable[ProtocolDAGResult]
            The :class:`.ProtocolDAGResult` objects to assemble aggregate
            quantities from.

        Returns
        -------
        ProtocolResult
            Aggregated results from many :class:`.ProtocolDAGResult` objects,
            all from a given :class:`.Protocol`.

        """
        return self.protocol.gather(protocol_dag_results=protocol_dag_results)

    def dump(self, file):
        """Dump this Transformation to a JSON file.

        Note that this is not space-efficient: for example, any ``Component``
        which is used in both ``ChemicalSystem`` objects will be represented
        twice in the JSON output.

        Parameters
        ----------
        file : Union[PathLike, FileLike]
            A pathlike of filelike to save this transformation to.
        """
        with ensure_filelike(file, mode="w") as f:
            json.dump(self.to_dict(), f, cls=JSON_HANDLER.encoder, sort_keys=True)

    @classmethod
    def load(cls, file):
        """Create a Transformation from a JSON file.

        Parameters
        ----------
        file : Union[PathLike, FileLike]
            A pathlike or filelike to read this transformation from.
        """
        with ensure_filelike(file, mode="r") as f:
            dct = json.load(f, cls=JSON_HANDLER.decoder)

        return cls.from_dict(dct)


class Transformation(TransformationBase):
    _stateA: ChemicalSystem
    _stateB: ChemicalSystem
    _name: Optional[str]
    _mapping: Optional[Union[ComponentMapping, list[ComponentMapping]]]
    _protocol: Protocol

    def __init__(
        self,
        stateA: ChemicalSystem,
        stateB: ChemicalSystem,
        protocol: Protocol,
        mapping: Optional[Union[ComponentMapping, list[ComponentMapping], dict[str, ComponentMapping]]] = None,
        name: Optional[str] = None,
    ):
        """Two chemical states with a method for estimating the free energy
        difference between them.

        Connects two :class:`.ChemicalSystem` objects, with directionality, and
        relates these to a :class:`.Protocol` which will provide an estimate of
        the free energy difference between these systems. Used as an edge of an
        :class:`.AlchemicalNetwork`.

        Parameters
        ----------
        stateA, stateB : ChemicalSystem
            The start (A) and end (B) states of the transformation.
        protocol : Protocol
            The method used to estimate the free energy difference between
            states A and B.
        mapping : Optional[Union[ComponentMapping, list[ComponentMapping]]]
            The details of any transformations between :class:`.Component` \s
            of the two states.
        name : str, optional
            A human-readable name for this transformation.

        """
        if isinstance(mapping, dict):
            warnings.warn(("mapping input as a dict is deprecated, "
                           "instead use either a single Mapping or list"),
                          DeprecationWarning)
            mapping = list(mapping.values())

        self._stateA = stateA
        self._stateB = stateB
        self._protocol = protocol
        self._mapping = mapping
        self._name = name

    def __repr__(self):
        return f"{self.__class__.__name__}(stateA={self.stateA}, stateB={self.stateB}, protocol={self.protocol}, name={self.name})"

    @property
    def stateA(self) -> ChemicalSystem:
        """The starting :class:`.ChemicalSystem` for the transformation."""
        return self._stateA

    @property
    def stateB(self) -> ChemicalSystem:
        """The ending :class:`.ChemicalSystem` for the transformation."""
        return self._stateB

    @property
    def protocol(self) -> Protocol:
        """The :class:`.Protocol` used to perform the transformation.

        This protocol estimates the free energy differences between ``stateA``
        and ``stateB`` :class:`.ChemicalSystem` objects. It includes all details
        needed to perform required simulations/calculations and encodes the
        alchemical pathway used.
        """
        return self._protocol

    @property
    def mapping(self) -> Optional[Union[ComponentMapping, list[ComponentMapping]]]:
        """The mappings relevant for this Transformation"""
        return self._mapping

    def _to_dict(self) -> dict:
        return {
            "stateA": self.stateA,
            "stateB": self.stateB,
            "protocol": self.protocol,
            "mapping": self.mapping,
            "name": self.name,
        }

    def create(
        self,
        *,
        extends: Optional[ProtocolDAGResult] = None,
        name: Optional[str] = None,
    ) -> ProtocolDAG:
        """
        Returns a ``ProtocolDAG`` executing this ``Transformation.protocol``.
        """
        return self.protocol.create(
            stateA=self.stateA,
            stateB=self.stateB,
            mapping=self.mapping,
            extends=extends,
            name=name,
            transformation_key=self.key,
        )


class NonTransformation(TransformationBase):

    def __init__(
        self,
        system: ChemicalSystem,
        protocol: Protocol,
        name: Optional[str] = None,
    ):
        """A non-alchemical edge of an alchemical network.

        A "transformation" that performs no transformation at all.
        Technically a self-loop, or an edge with the same ``ChemicalSystem`` at
        either end.

        Functionally used for applying a dynamics protocol to a ``ChemicalSystem``
        that performs no alchemical transformation at all. This allows e.g.
        equilibrium MD to be performed on a ``ChemicalSystem`` as desired alongside
        alchemical protocols between it and and other ``ChemicalSystem`` objects.

        Parameters
        ----------
        system : ChemicalSystem
            The (identical) end states of the "transformation" to be sampled
        protocol : Protocol
            The sampling method to use on the ``system``
        name : str, optional
            A human-readable name for this transformation.

        """

        self._system = system
        self._protocol = protocol
        self._name = name

    def __repr__(self):
        return f"{self.__class__.__name__}(system={self.system}, protocol={self.protocol}, name={self.name})"

    @property
    def stateA(self) -> ChemicalSystem:
        """The :class:`.ChemicalSystem` this "transformation" samples.

        Synonomous with ``system`` attribute.

        """
        return self._system

    @property
    def stateB(self) -> ChemicalSystem:
        """The :class:`.ChemicalSystem` this "transformation" samples.

        Synonomous with ``system`` attribute.

        """
        return self._system

    @property
    def system(self) -> ChemicalSystem:
        """The :class:`.ChemicalSystem` this "transformation" samples."""
        return self._system

    @property
    def protocol(self):
        """The :class:`.Protocol` for sampling dynamics of the ``system``.

        Includes all details needed to perform required
        simulations/calculations.
        """
        return self._protocol

    def _to_dict(self) -> dict:
        return {
            "system": self.system,
            "protocol": self.protocol,
            "name": self.name,
        }

    def create(
        self,
        *,
        extends: Optional[ProtocolDAGResult] = None,
        name: Optional[str] = None,
    ) -> ProtocolDAG:
        """
        Returns a ``ProtocolDAG`` executing this ``NonTransformation.protocol``.
        """
        # TODO: once we have an implicit component mapping concept, use this
        # here instead of None to allow use of alchemical protocols with
        # NonTransformations
        return self.protocol.create(
            stateA=self.system,
            stateB=self.system,
            mapping=None,
            extends=extends,
            name=name,
            transformation_key=self.key,
        )
