# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from .protocol import Protocol


class Transformation:
    """An edge of an alchemical network.

    Connects two `ChemicalState`s, with directionality.

    Attributes
    ----------
    start : ChemicalState
    end : ChemicalState
    protocol : Protocol
        The protocol used to perform the transformation.
        Includes all details needed to perform required
        simulations/calculations and encodes the alchemical pathway used.

    """

    def __init__(
            self,
            start,
            end,
            protocol
            ):

        self._protocol = protocol

    def start(self):
        return self._start

    def end(self):
        return self._end


    @property
    def protocol(self):
        """The protocol for sampling the transformation to derive free energy
        differences between the `ChemicalState`s on either end.

        """
        return self._protocol
