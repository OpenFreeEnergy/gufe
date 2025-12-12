# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import abc
from collections.abc import Iterator

import gufe

from ..tokenization import GufeTokenizable
from . import LigandAtomMapping
from .atom_mapping import AtomMapping


class AtomMapper(GufeTokenizable):
    """A class for manufacturing mappings

    Implementations of this class can require an arbitrary and non-standardised
    number of input arguments to create.

    Implementations of this class provide the :meth:`.suggest_mappings` method
    """

    @abc.abstractmethod
    def suggest_mappings(self, A: gufe.Component, B: gufe.Component) -> Iterator[AtomMapping]:
        """Suggests possible mappings between two Components

        Suggests zero or more :class:`.AtomMapping` objects, which are possible
        atom mappings between two :class:`.Component` objects.
        """
        ...


class LigandAtomMapper(gufe.AtomMapper):
    """
    Suggest atom mappings between two :class:`SmallMoleculeComponent` instances.

    Subclasses will typically implement the ``_mappings_generator`` method,
    which returns an iterable of :class:`.LigandAtomMapping` suggestions.
    """

    @abc.abstractmethod
    def _mappings_generator(
        self,
        componentA: SmallMoleculeComponent,
        componentB: SmallMoleculeComponent,
    ) -> Iterable[dict[int, int]]:
        """
        Suggest mapping options for the input molecules.

        Parameters
        ----------
        componentA, componentB : rdkit.Mol
            the two molecules to create a mapping for

        Returns
        -------
        Iterable[dict[int, int]] :
            an iterable over proposed mappings from componentA to componentB
        """
        ...

    def suggest_mappings(
        self,
        componentA: SmallMoleculeComponent,
        componentB: SmallMoleculeComponent,
    ) -> Iterable[LigandAtomMapping]:
        """
        Suggest :class:`.LigandAtomMapping` options for the input molecules.

        Parameters
        ---------
        componentA, componentB : :class:`.SmallMoleculeComponent`
            the two molecules to create a mapping for

        Returns
        -------
        Iterable[LigandAtomMapping] :
            an iterable over proposed mappings
        """
        # For this base class, implementation is redundant with
        # _mappings_generator. However, we keep it separate so that abstract
        # subclasses of this can customize suggest_mappings while always
        # maintaining the consistency that concrete implementations must
        # implement _mappings_generator.

        for map_dct in self._mappings_generator(componentA, componentB):
            yield LigandAtomMapping(componentA, componentB, map_dct)
