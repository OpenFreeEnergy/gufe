# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import abc
from collections.abc import Iterator
import gufe


class AtomMapper(abc.ABC):
    """A class for manufacturing mappings"""
    def suggest_mappings(self, A: gufe.Component,
                         B: gufe.Component) -> Iterator[gufe.AtomMapping]:
        ...
