# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from typing import Union

from ..base import GufeTokenizable


class Component(GufeTokenizable):
    """Base class for members of a ChemicalSystem"""

    def __repr__(self):
        return f"{self.__class__.__name__}(name={self.name})"

    @property
    @abc.abstractmethod
    def name(self) -> str:
        pass

    @abc.abstractmethod
    def _to_dict(self) -> dict:
        ...

    @classmethod
    @abc.abstractmethod
    def _from_dict(cls, d: dict):
        ...

    @property
    @abc.abstractmethod
    def total_charge(self) -> Union[int, None]:
        """Net formal charge for the Component if defined"""
        ...
