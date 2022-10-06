# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from collections.abc import Mapping, Iterable


import gufe
from .componentmapping import ComponentMapping


class AtomMapping(ComponentMapping, abc.ABC):
    _molA: gufe.Component
    _molB: gufe.Component

    """A mapping between two different atom-based Components"""
    @abc.abstractmethod
    def _to_dict(self) -> dict:
        ...

    @classmethod
    @abc.abstractmethod
    def _from_dict(cls, d: dict):
        ...

    @property
    @abc.abstractmethod
    def molA(self) -> gufe.Component:
        """A copy of the first Component in the mapping"""
        return self._molA

    @property
    @abc.abstractmethod
    def molB(self) -> gufe.Component:
        """A copy of the second Component in the mapping"""
        return self._molB

    @property
    @abc.abstractmethod
    def molA_to_molB(self) -> Mapping[int, int]:
        """Maps the index of an item from Component A onto Component B

        Not all indices will be resolvable, these items are not mapped,
        resulting in a KeyError
        """
        ...

    @property
    @abc.abstractmethod
    def molB_to_molA(self) -> Mapping[int, int]:
        ...

    @property
    @abc.abstractmethod
    def molA_unique(self) -> Iterable[int]:
        """Indices of atoms in mol A that aren't mappable to B"""
        ...

    @property
    @abc.abstractmethod
    def molB_unique(self) -> Iterable[int]:
        """Indices of atoms in mol B that aren't mappable to A"""
        ...
