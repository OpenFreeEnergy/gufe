# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from typing import Optional

from openff.toolkit.utils.serialization import Serializable

from .chemicalsystem import ChemicalSystem
from .protocols import Protocol
from .mapping import AtomMapping


class Transformation(Serializable):
    """An edge of an alchemical network.

    Connects two `ChemicalSystem`s, with directionality.

    Attributes
    ----------
    initial : ChemicalSystem
    final : ChemicalSystem
    protocol : Protocol
        The protocol used to perform the transformation.
        Includes all details needed to perform required
        simulations/calculations and encodes the alchemical pathway used.
    identifier
        Optional identifier for the transformation; set this to a unique value
        if adding multiple, otherwise identical transformations to the same
        `AlchemicalNetwork` to avoid deduplication

    """

    def __init__(
            self, 
            initial: ChemicalSystem,
            final: ChemicalSystem,
            mapping: Optional[AtomMapping] = None,
            protocol: Optional[Protocol] = None,
            identifier: Optional[str] = None,
        ):

        self._initial = initial
        self._final = final
        self._mapping = mapping
        self._identifier = identifier

        self._protocol = protocol

    def __repr__(self):
        return f"{self.__class__.__name__}(initial={self.initial}, final={self.final}, protocol={self.protocol})"

    @property
    def initial(self):
        return self._initial

    @property
    def final(self):
        return self._final

    @property
    def protocol(self):
        """The protocol for sampling the transformation to derive free energy
        differences between `initial` and `final` `ChemicalSystem`s.

        """
        return self._protocol

    @property
    def mapping(self):
        """The mapping between atoms in `initial` to `final` `ChemicalSystem`s for
        the transformation.

        """
        return self._mapping

    @property
    def identifier(self):
        return self._identifier

    def __lt__(self, other):
        return hash(self) < hash(other)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        if self._identifier != other.identifier:
            return False
        if self._mapping != other.mapping:
            return False
        if self._final != other.final:
            return False
        if self._initial != other.initial:
            return False

        return True

    def __hash__(self):
        return hash(
            (
                self._initial,
                self._final,
                self._mapping,
                self._identifier,
            )
        )

    # TODO: broken without a `Protocol` registry of some kind
    # should then be changed to use the registry
    def to_dict(self) -> dict:
        return {
            "initial": self.initial.to_dict(),
            "final": self.final.to_dict(),
            "protocol": self.protocol.to_dict(),
        }

    # TODO: broken without a `Protocol` registry of some kind
    # should then be changed to use the registry
    @classmethod
    def from_dict(cls, d: dict):
        return cls(
                initial=ChemicalSystem.from_dict(d['initial']),
                final=ChemicalSystem.from_dict(d['final']),
                protocol=Protocol.from_dict(d['protocol'])
                )

    def run(self):
        return self.protocol.run(
                                 initial=self.initial, 
                                 final=self.final,
                                 mapping=self.mapping
                                )


# we subclass `Transformation` here for typing simplicity
class NonTransformation(Transformation):
    """A non-alchemical edge of an alchemical network.

    A "transformation" that performs no transformation at all.
    Technically a self-loop, or an edge with the same `ChemicalSystem` at either final.

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

    """

    def __init__(
            self,
            chemicalsystem: ChemicalSystem,
            protocol: Optional[Protocol] = None
        ):

        self._chemicalsystem = chemicalsystem
        self._protocol = protocol

    @property
    def initial(self):
        return self._chemicalsystem

    @property
    def final(self):
        return self._chemicalsystem

    @property
    def system(self):
        return self._chemicalsystem

    @property
    def protocol(self):
        """The protocol for sampling dynamics of the `ChemicalSystem`."""
        return self._protocol

    # TODO: broken without a `Protocol` registry of some kind
    # should then be changed to use the registry
    def to_dict(self) -> dict:
        return {
            "chemicalsystem": self.chemicalsystem.to_dict(),
            "protocol": self.protocol.to_dict(),
        }

    # TODO: broken without a `Protocol` registry of some kind
    # should then be changed to use the registry
    @classmethod
    def from_dict(cls, d: dict):
        return cls(
                chemicalsystem=ChemicalSystem.from_dict(d['chemicalsystem']),
                protocol=Protocol.from_dict(d['protocol'])
                )
