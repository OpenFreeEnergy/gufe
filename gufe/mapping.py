# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from openff.toolkit.utils.serialization import Serializable


class Mapping(Serializable):
    """A mapping, usually of atoms, between two `ChemicalSystem`s.

    Some, but not all, transformation protocols require a mapping to execute.

    """

    def to_dict(self) -> dict:
        ...

    @classmethod
    def from_dict(cls, d: dict):
        ...
