# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import abc
from collections.abc import Iterator
import gufe

from .atom_mapping import AtomMapping


class AtomMapper(abc.ABC):
    """A class for manufacturing mappings"""
    @abc.abstractmethod
    def suggest_mappings(self, A: gufe.Component,
                         B: gufe.Component) -> Iterator[AtomMapping]:
        ...
