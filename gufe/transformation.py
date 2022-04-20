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
    start : ChemicalSystem
    end : ChemicalSystem
    protocol : Protocol
        The protocol used to perform the transformation.
        Includes all details needed to perform required
        simulations/calculations and encodes the alchemical pathway used.

    """

    def __init__(
            self, 
            start: ChemicalSystem,
            end: ChemicalSystem,
            mapping: Optional[AtomMapping] = None,
            protocol: Optional[Protocol] = None
        ):

        self._start = start
        self._end = end
        self._mapping = mapping

        self._protocol = protocol

    def __repr__(self):
        return f"{self.__class__.__name__}(start={self.start}, end={self.end}, protocol={self.protocol})"

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def protocol(self):
        """The protocol for sampling the transformation to derive free energy
        differences between the `ChemicalSystem`s on either end.

        """
        return self._protocol

    @property
    def mapping(self):
        """The mapping 

        """
        return self._protocol

    def __hash__(self):
        return hash(
            (
                tuple(sorted(self._components.items())),
                self._box_vectors.tobytes(),
                self._identifier,
            )
        )

    # TODO: broken without a `Protocol` registry of some kind
    # should then be changed to use the registry
    def to_dict(self) -> dict:
        return {
            "start": self.start.to_dict(),
            "end": self.end.to_dict(),
            "protocol": self.protocol.to_dict(),
        }

    # TODO: broken without a `Protocol` registry of some kind
    # should then be changed to use the registry
    @classmethod
    def from_dict(cls, d: dict):
        return cls(
                start=ChemicalSystem.from_dict(d['start']),
                end=ChemicalSystem.from_dict(d['end']),
                protocol=Protocol.from_dict(d['protocol'])
                )

    def run(self):
        return self.protocol.run(
                                 start=self.start, 
                                 end=self.end,
                                 mapping=self.mapping
                                )

# we subclass `Transformation` here for typing simplicity
class NonTransformation(Transformation):
    """A non-alchemical edge of an alchemical network.

    A "transformation" that performs no transformation at all.
    Technically a self-loop, or an edge with the same `ChemicalSystem` at either end.

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
    def start(self):
        return self._chemicalsystem

    @property
    def end(self):
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
