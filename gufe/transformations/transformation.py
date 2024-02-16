# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import json
from collections.abc import Iterable
from typing import Optional, Union

from ..chemicalsystem import ChemicalSystem
from ..mapping import ComponentMapping
from ..protocols import Protocol, ProtocolDAG, ProtocolDAGResult, ProtocolResult
from ..tokenization import JSON_HANDLER, GufeTokenizable
from ..utils import ensure_filelike


class Transformation(GufeTokenizable):
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
        mapping: Optional[Union[ComponentMapping, list[ComponentMapping]]] = None,
        name: Optional[str] = None,
    ):
        r"""Two chemical states with a method for estimating free energy difference

        Connects two :class:`.ChemicalSystem` objects, with directionality,
        and relates this to a :class:`.Protocol` which will provide an estimate of
        the free energy difference of moving between these systems.
        Used as an edge of an :class:`.AlchemicalNetwork`.

        Parameters
        ----------
        stateA, stateB: ChemicalSystem
           The start (A) and end (B) states of the transformation
        protocol: Protocol
           The method used to estimate the free energy difference between states
           A and B
        mapping : Optional[Union[ComponentMapping, list[ComponentMapping]]]
           the details of any transformations between :class:`.Component` \s of
           the two states
        name : str, optional
           a human-readable tag for this transformation
        """
        self._stateA = stateA
        self._stateB = stateB
        self._mapping = mapping
        self._name = name

        self._protocol = protocol

    @classmethod
    def _defaults(cls):
        return super()._defaults()

    def __repr__(self):
        return f"{self.__class__.__name__}(stateA={self.stateA}, " f"stateB={self.stateB}, protocol={self.protocol})"

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
        """The protocol used to perform the transformation.

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

    @property
    def name(self) -> Optional[str]:
        """
        Optional identifier for the transformation; used as part of its hash.

        Set this to a unique value if adding multiple, otherwise identical
        transformations to the same :class:`AlchemicalNetwork` to avoid
        deduplication.
        """
        return self._name

    def _to_dict(self) -> dict:
        return {
            "stateA": self.stateA,
            "stateB": self.stateB,
            "protocol": self.protocol,
            "mapping": self.mapping,
            "name": self.name,
        }

    @classmethod
    def _from_dict(cls, d: dict):
        return cls(**d)

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

    def gather(self, protocol_dag_results: Iterable[ProtocolDAGResult]) -> ProtocolResult:
        """
        Gather multiple ``ProtocolDAGResult`` into a single ``ProtocolResult``.

        Parameters
        ----------
        protocol_dag_results : Iterable[ProtocolDAGResult]
            The ``ProtocolDAGResult`` objects to assemble aggregate quantities
            from.

        Returns
        -------
        ProtocolResult
            Aggregated results from many ``ProtocolDAGResult`` objects, all from
            a given ``Protocol``.

        """
        return self.protocol.gather(protocol_dag_results=protocol_dag_results)

    def dump(self, file):
        """Dump this Transformation to a JSON file.

        Note that this is not space-efficient: for example, any
        ``Component`` which is used in both ``ChemicalSystem`` objects will be
        represented twice in the JSON output.

        Parameters
        ----------
        file : Union[PathLike, FileLike]
            a pathlike of filelike to save this transformation to.
        """
        with ensure_filelike(file, mode="w") as f:
            json.dump(self.to_dict(), f, cls=JSON_HANDLER.encoder, sort_keys=True)

    @classmethod
    def load(cls, file):
        """Create a Transformation from a JSON file.

        Parameters
        ----------
        file : Union[PathLike, FileLike]
            a pathlike or filelike to read this transformation from
        """
        with ensure_filelike(file, mode="r") as f:
            dct = json.load(f, cls=JSON_HANDLER.decoder)

        return cls.from_dict(dct)


# we subclass `Transformation` here for typing simplicity
class NonTransformation(Transformation):
    """A non-alchemical edge of an alchemical network.

    A "transformation" that performs no transformation at all.
    Technically a self-loop, or an edge with the same ``ChemicalSystem`` at
    either end.

    Functionally used for applying a dynamics protocol to a ``ChemicalSystem``
    that performs no alchemical transformation at all. This allows e.g.
    equilibrium MD to be performed on a ``ChemicalSystem`` as desired alongside
    alchemical protocols between it and and other ``ChemicalSystem`` objects.
    """

    def __init__(
        self,
        system: ChemicalSystem,
        protocol: Protocol,
        name: Optional[str] = None,
    ):

        self._system = system
        self._name = name
        self._protocol = protocol

    @property
    def stateA(self):
        return self._system

    @property
    def stateB(self):
        return self._system

    @property
    def system(self) -> ChemicalSystem:
        return self._system

    @property
    def protocol(self):
        """
        The protocol for sampling dynamics of the `ChemicalSystem`.

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

    @classmethod
    def _from_dict(cls, d: dict):
        return cls(**d)

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
            stateA=self.system,
            stateB=self.system,
            mapping=None,
            extends=extends,
            name=name,
            transformation_key=self.key,
        )
