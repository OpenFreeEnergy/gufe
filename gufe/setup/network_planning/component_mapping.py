# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import abc

import gufe
from gufe.tokenization import GufeTokenizable


class ComponentMapping(GufeTokenizable, abc.ABC):
    """A relationship between two Components stating that they transform in some way

    For components that are atom-based is specialised to :class:`.AtomMapping`
    """
    _componentA: gufe.Component
    _componentB: gufe.Component

    def __init__(self, componentA: gufe.Component, componentB: gufe.Component):
        self._componentA = componentA
        self._componentB = componentB
        # self.componentA_to_componentB # TODO: is that something we want here, thinking beyond AtomMappings?

    def __contains__(self, item: gufe.Component):
        return item == self._componentA or item == self._componentB

    @property
    def componentA(self) -> gufe.Component:
        """A copy of the first Component in the mapping"""
        return self._componentA

    @property
    def componentB(self) -> gufe.Component:
        """A copy of the second Component in the mapping"""
        return self._componentB