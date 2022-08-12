# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from typing import Optional, Iterable

from openff.toolkit.utils.serialization import Serializable
from ..tokenize import GufeTokenizable

from ..chemicalsystem import ChemicalSystem
from ..protocols import Protocol, ProtocolDAG, ProtocolResult, ProtocolDAGResult
from ..mapping import Mapping


class Transformation:  #(GufeTokenizable):
    """An edge of an alchemical network.

    Connects two `ChemicalSystem`s, with directionality.

    Attributes
    ----------
    stateA : ChemicalSystem
        The starting `ChemicalSystem` for the transformation.
    stateB : ChemicalSystem
        The ending `ChemicalSystem` for the transformation.
    protocol : Protocol
        The protocol used to perform the transformation.
        Includes all details needed to perform required
        simulations/calculations and encodes the alchemical pathway used.
        May also correspond to an experimental result.
    mapping : Optional[Mapping]
        Mapping of e.g. atoms between the `stateA` and `stateB`
        `ChemicalSystem`s.
    name : Optional[str]
        Optional identifier for the transformation; set this to a unique value
        if adding multiple, otherwise identical transformations to the same
        `AlchemicalNetwork` to avoid deduplication

    """

    def __init__(
        self,
        stateA: ChemicalSystem,
        stateB: ChemicalSystem,
        protocol: Protocol,
        mapping: Optional[Mapping] = None,
        name: Optional[str] = None,
    ):

        self._stateA = stateA
        self._stateB = stateB
        self._mapping = mapping
        self._name = name

        self._protocol = protocol

    def _defaults(self):
        return super().defaults()

    def __repr__(self):
        return f"{self.__class__.__name__}(stateA={self.stateA}, "\
               f"stateB={self.stateB}, protocol={self.protocol})"

    @property
    def stateA(self):
        """The starting `ChemicalSystem` for the transformation."""
        return self._stateA

    @property
    def stateB(self):
        """The ending `ChemicalSystem` for the transformation."""
        return self._stateB

    @property
    def protocol(self):
        """The protocol for sampling the transformation to derive free energy
        differences between `stateA` and `stateB` `ChemicalSystem`s.

        """
        return self._protocol

    @property
    def mapping(self):
        """The mapping between atoms in `stateA` to `stateB`"""
        return self._mapping

    @property
    def name(self):
        """User-specified for the transformation; used as part of its hash."""
        return self._name

    def __lt__(self, other):
        return hash(self) < hash(other)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        if self._name != other.name:
            return False
        if self._protocol != other.protocol:
            return False
        if self._mapping != other.mapping:
            return False
        if self._stateB != other.stateB:
            return False
        if self._stateA != other.stateA:
            return False

        return True

    def __hash__(self):
        return hash(
            (
                self._stateA,
                self._stateB,
                self._mapping,
                self._name,
                self._protocol,
            )
        )

    # TODO: broken without a `Protocol` registry of some kind
    # should then be changed to use the registry
    def _to_dict(self) -> dict:
        return {
            "stateA": self.stateA.to_dict(),
            "stateB": self.stateB.to_dict(),
            "protocol": self.protocol.to_dict(),
            "mapping": self.mapping.to_dict(),
        }

    # TODO: broken without a `Protocol` registry of some kind
    # should then be changed to use the registry
    @classmethod
    def _from_dict(cls, d: dict):
        return cls(
            stateA=ChemicalSystem.from_dict(d["stateA"]),
            stateB=ChemicalSystem.from_dict(d["stateB"]),
            protocol=Protocol.from_dict(d["protocol"]),
            mapping=Mapping.from_dict(d["mapping"]),
        )

    def create(self) -> ProtocolDAG:
        """Returns a `ProtocolDAG` executing this `Transformation.protocol`."""
        return self.protocol.create(
            stateA=self.stateA,
            stateB=self.stateB,
            mapping=self.mapping,
            name=str(self.__hash__()),
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
        return self.protocol.gather(protocol_dag_results=protocol_dag_results)


# we subclass `Transformation` here for typing simplicity
class NonTransformation(Transformation):
    """A non-alchemical edge of an alchemical network.

    A "transformation" that performs no transformation at all.
    Technically a self-loop, or an edge with the same `ChemicalSystem` at
    either end.

    Functionally used for applying a dynamics protocol to a `ChemicalSystem`
    that performs no alchemical transformation at all. This allows e.g.
    equilibrium MD to be performed on a `ChemicalSystem` as desired alongside
    alchemical protocols between it and and other `ChemicalSystem`s.

    Attributes
    ----------
    system : ChemicalSystem
    protocol : Protocol
        The protocol used to perform the dynamics. Includes all details needed
        to perform required simulations/calculations.
    name : Optional[str]
        Optional identifier for the nontransformation; set this to a unique
        value if adding multiple, otherwise identical transformations to the
        same `AlchemicalNetwork` to avoid deduplication

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
    def system(self):
        return self._system

    @property
    def protocol(self):
        """The protocol for sampling dynamics of the `ChemicalSystem`."""
        return self._protocol

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        if self._name != other.name:
            return False
        if self._protocol != other.protocol:
            return False
        if self._system != other.system:
            return False

        return True

    def __hash__(self):
        return hash(
            (
                self._name,
                self._protocol,
                self._system,
            )
        )

    # TODO: broken without a `Protocol` registry of some kind
    # should then be changed to use the registry
    def _to_dict(self) -> dict:
        return {
            "system": self.system.to_dict(),
            "protocol": self.protocol.to_dict(),  # TODO: include defaults?
        }

    # TODO: broken without a `Protocol` registry of some kind
    # should then be changed to use the registry
    @classmethod
    def _from_dict(cls, d: dict):
        return cls(
            system=ChemicalSystem.from_dict(d["system"]),
            protocol=Protocol.from_dict(d["protocol"]),
        )

    def create(self) -> ProtocolDAG:
        """Returns a `ProtocolDAG` executing this `Transformation.protocol`."""
        return self.protocol.create(
            stateA=self.system, stateB=self.system, name=str(self.__hash__())
        )
