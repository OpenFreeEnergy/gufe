# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from typing import Union

from ..tokenization import GufeTokenizable


class Component(GufeTokenizable):
    """Base class for an individual item in a larger chemical system

    Individual items could be as small as a single molecule, or as large as an
    entire biological assembly.

    These items are used as pieces of a :class:`ChemicalSystem`.
    """

    def __repr__(self):
        return f"{self.__class__.__name__}(name={self.name})"

    @property
    @abc.abstractmethod
    def name(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def total_charge(self) -> Union[int, None]:
        """Net formal charge for the ``Component``, if defined."""
        ...
