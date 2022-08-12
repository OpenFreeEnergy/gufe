# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import abc
from typing import Iterable

from gufe import ChemicalSystem, Component
from gufe.base import GufeTokenizable
from . import ComponentMapping


class SystemMapping(GufeTokenizable, abc.ABC):
    """Class defining mappings between two Systems

    Often contains many mappings
    """
    stateA: ChemicalSystem
    stateB: ChemicalSystem
    mappings: list[ComponentMapping]

    def __init__(self, stateA: ChemicalSystem, stateB: ChemicalSystem,
                 mappings: Iterable[ComponentMapping]):
        self._stateA = stateA
        self._stateB = stateB
        self._mappings = tuple(mappings)
        # TODO: what sanity checks here are needed?

    def is_mapped(self, thing: Component) -> bool:
        """Does a given Component have an associated mapping"""
        return any(thing in m for m in self._mappings)

    def get_mapping(self, thing: Component) -> ComponentMapping:
        """Return the mapping in which a Component features

        Raises
        ------
        ValueError
          if no mapping found
        """
        for m in self._mappings:
            if thing in m:
                return m
        else:
            raise ValueError("No mapping found!")

