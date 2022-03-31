from .protocols import Protocol


class Transformation:
    """An edge of an alchemical network.

    Connects two chemical states, with directionality.

    Attributes
    ----------
    protocol : Protocol
        The protocol used to perform the transformation.
        Includes all details needed to perform calculation and encodes the
        alchemical pathway used.

    """

    def __init__(
            self,
            protocol
            ):

        self._protocol = protocol


    # these are ambiguous; will need a better way to consider this
    # these could be expensive operations
    def dg(self, estimator=None):
        """Free energy difference of transformation based on given estimator,
        using only data available for this edge.

        """
        ...

    def ddg(self, estimator=None):
        ...

    def chemicalstate_start(self):
        ...

    def chemicalstate_end(self):
        ...


    @property
    def protocol(self):
        """The protocol for sampling the transformation to derive free energy
        differences between the `ChemicalState`s on either end.

        """
        return self._protocol
