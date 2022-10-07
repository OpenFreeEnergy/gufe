# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import abc

import gufe
from gufe.tokenization import GufeTokenizable


class ComponentMapping(GufeTokenizable, abc.ABC):
    _molA: gufe.Component
    _molB: gufe.Component

    def __contains__(self, item: gufe.Component):
        return item == self._molA or item == self._molB
