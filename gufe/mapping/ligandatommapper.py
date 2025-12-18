import abc
from typing import Iterator

from gufe import SmallMoleculeComponent

from .atom_mapper import AtomMapper
from .ligandatommapping import LigandAtomMapping


class LigandAtomMapper(AtomMapper):
    """
    Suggest atom mappings between two :class:`SmallMoleculeComponent` instances.

    Subclasses will typically implement the ``_mappings_generator`` method,
    which returns an Iterator of :class:`.LigandAtomMapping` suggestions.
    """

    @abc.abstractmethod
    def _mappings_generator(
        self,
        componentA: SmallMoleculeComponent,
        componentB: SmallMoleculeComponent,
    ) -> Iterator[dict[int, int]]:
        """
        Suggest mapping options for the input molecules.

        Parameters
        ----------
        componentA, componentB : rdkit.Mol
            the two molecules to create a mapping for

        Returns
        -------
        Iterator[dict[int, int]] :
            an Iterator over proposed mappings from componentA to componentB
        """
        ...

    def suggest_mappings(
        self,
        # TODO: fix overrides when we move to min python 3.12 - see https://peps.python.org/pep-0695/#summary-examples
        componentA: SmallMoleculeComponent,  # type: ignore[override]
        componentB: SmallMoleculeComponent,  # type: ignore[override]
    ) -> Iterator[LigandAtomMapping]:
        """
        Suggest :class:`.LigandAtomMapping` options for the input molecules.

        Parameters
        ---------
        componentA, componentB : :class:`.SmallMoleculeComponent`
            the two molecules to create a mapping for

        Returns
        -------
        Iterator[LigandAtomMapping] :
            an Iterator over proposed mappings
        """
        # For this base class, implementation is redundant with
        # _mappings_generator. However, we keep it separate so that abstract
        # subclasses of this can customize suggest_mappings while always
        # maintaining the consistency that concrete implementations must
        # implement _mappings_generator.

        for map_dct in self._mappings_generator(componentA, componentB):
            yield LigandAtomMapping(componentA, componentB, map_dct)
