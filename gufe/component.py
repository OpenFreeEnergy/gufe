# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from typing import Union

from .base import GufeTokenizableMixin


COMPONENT_REGISTRY = {}

def register_component(cls):
    COMPONENT_REGISTRY[cls.__name__] = cls

    
class _ComponentMeta(type):
    def __new__(cls, name, bases, classdict):
        componentcls = super().__new__(cls, name, bases, classdict)
        register_component(componentcls)
        return componentcls


class _ABCComponentMeta(_ComponentMeta, abc.ABCMeta):
    # required to make use of abc.ABC in classes that use _ComponentMeta metaclass
    # see https://stackoverflow.com/questions/57349105/python-abc-inheritance-with-specified-metaclass 
    ...
        

class Component(abc.ABC, GufeTokenizableMixin, metaclass=_ABCComponentMeta):
    """Base class for members of a ChemicalSystem"""

    def __init__(self):

        # since __hash__ isn't inheritable, monkey-patching is the best we can
        # do besides defining it in every subclass
        self.__class__.__hash__ = Component.__hash__

    def __repr__(self):
        return f"{self.__class__.__name__}(name={self.name})"

    def __lt__(self, other):
        return hash(self) < hash(other)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented("Component comparisons must be of the same type")

        d = self.to_dict()
        do = other.to_dict()

        if set(d.keys()) == set(do.keys()):
            for key in d:
                if d[key] != do[key]:
                    return False

            return True

    def __hash__(self):
        return hash(self.key)

    def _gufe_tokenize(self):
        """Return a list of normalized inputs for `gufe.base.tokenize`.

        """
        return sorted(self.to_dict(include_defaults=False).items(), key=str)

    @property
    @abc.abstractmethod
    def name(self) -> str:
        pass

    def to_dict(self, include_defaults=True, keyencode_dependencies=False) -> dict:
        """Serialize to dict representation"""

        d = self._to_dict()
        d['__class__'] = self.__class__.__name__

        if keyencode_dependencies:
            d = self._keyencode_dependencies(d)

        if not include_defaults:
            for key, value in self.defaults.items():
                if d.get(key) == value:
                    d.pop(key)

        return d

    @abc.abstractmethod
    def _to_dict(self) -> dict:
        ...


    def from_dict(cls, d: dict, keydecode_dependencies=False):
        """Deserialize from dict representation"""
        component_type = d.pop('__class__', None)
        if component_type is None:
            raise KeyError("No `__class__` key found; unable to deserialize Component")

        return COMPONENT_REGISTRY[component_type]._from_dict(d)

    @abc.abstractmethod
    def _from_dict(cls, d: dict):
        ...

    @property
    @abc.abstractmethod
    def total_charge(self) -> Union[int, None]:
        """Net formal charge for the Component if defined"""
        ...
