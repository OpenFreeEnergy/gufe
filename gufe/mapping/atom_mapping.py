# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from collections.abc import Mapping, Iterable


import gufe
from .componentmapping import ComponentMapping


class AtomMapping(ComponentMapping, abc.ABC):
    _componentA: gufe.Component
    _componentB: gufe.Component

    """A mapping between two different atom-based Components"""

    @property
    def componentA(self) -> gufe.Component:
        """A copy of the first Component in the mapping"""
        return self._componentA

    @property
    def componentB(self) -> gufe.Component:
        """A copy of the second Component in the mapping"""
        return self._componentB

    @property
    @abc.abstractmethod
    def componentA_to_componentB(self) -> Mapping[int, int]:
        """Maps the index of an item from Component A onto Component B

        Not all indices will be resolvable, these items have no corresponding
        entity in the other component (e.g. the atom disappears), therefore
        resulting in a KeyError on query
        """
        ...

    @property
    @abc.abstractmethod
    def componentB_to_componentA(self) -> Mapping[int, int]:
        """Similar to A to B, but reversed."""
        ...

    @property
    @abc.abstractmethod
    def componentA_unique(self) -> Iterable[int]:
        """Indices of atoms in component A that aren't mappable to B"""
        ...

    @property
    @abc.abstractmethod
    def componentB_unique(self) -> Iterable[int]:
        """Indices of atoms in component B that aren't mappable to A"""
        ...
