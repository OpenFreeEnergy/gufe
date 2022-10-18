# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import abc

import gufe
from gufe.tokenization import GufeTokenizable


class ComponentMapping(GufeTokenizable, abc.ABC):
    _componentA: gufe.Component
    _componentB: gufe.Component

    def __contains__(self, item: gufe.Component):
        return item == self._componentA or item == self._componentB
