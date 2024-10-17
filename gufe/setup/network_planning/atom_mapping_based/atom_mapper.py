# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import abc
from collections.abc import Iterator
import gufe

from ..component_mapper import ComponentMapper
from .atom_mapping import AtomMapping


class AtomMapper(ComponentMapper):
    """A class for manufacturing mappings

    Implementations of this class can require an arbitrary and non-standardised
    number of input arguments to create.

    Implementations of this class provide the :meth:`.suggest_mappings` method
    """

    @abc.abstractmethod
    def suggest_mappings(self,
                         A: gufe.Component,
                         B: gufe.Component
                         ) -> Iterator[AtomMapping]:
        """Suggests possible mappings between two Components

        Suggests zero or more :class:`.AtomMapping` objects, which are possible
        atom mappings between two :class:`.Component` objects.
        """
        raise NotImplementedError("This function was not implemented.")

