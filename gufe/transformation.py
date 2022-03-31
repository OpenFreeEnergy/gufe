# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from .protocol import Protocol


class Transformation:
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

    def __init__(self, start, end, protocol):

        self._start = start
        self._end = end

        self._protocol = protocol

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


class NonTransformation:
    """A non-alchemical edge of an alchemical network.

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

    def __init__(self, chemicalsystem, protocol):

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
