# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from collections.abc import Mapping
from openff.toolkit.utils.serialization import Serializable

import gufe


class AtomMapping(Serializable, abc.ABC):
    """A mapping between two different Components"""
    @abc.abstractmethod
    def to_dict(self) -> dict:
        ...

    @classmethod
    @abc.abstractmethod
    def from_dict(cls, d: dict):
        ...

    @property
    @abc.abstractmethod
    def molA(self) -> gufe.Component:
        """A reference to the first Component in the mapping"""
        ...

    @property
    @abc.abstractmethod
    def molB(self) -> gufe.Component:
        """A reference to the second Component in the mapping"""
        ...

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
