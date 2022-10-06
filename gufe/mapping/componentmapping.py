# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import abc

import gufe
from gufe.tokenization import GufeTokenizable


class ComponentMapping(GufeTokenizable, abc.ABC):
    stateA: gufe.Component
    stateB: gufe.Component

    def __contains__(self, item: gufe.Component):
        return item == self.stateA or item == self.stateB
