# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from typing import TypeVar

import gufe
from gufe.tokenization import GufeTokenizable

T = TypeVar("T", bound=gufe.Component)


class ComponentMapping(GufeTokenizable):
    """A relationship between two Components stating that they transform in some way

    For components that are atom-based is specialised to :class:`.AtomMapping`
    """

    _componentA: T
    _componentB: T

    def __init__(self, componentA: T, componentB: T):
        self._componentA = componentA
        self._componentB = componentB

    def __contains__(self, item: T):
        return item == self._componentA or item == self._componentB

    @property
    def componentA(self) -> T:
        """The first Component in the mapping"""
        return self._componentA

    @property
    def componentB(self) -> T:
        """The second Component in the mapping"""
        return self._componentB
