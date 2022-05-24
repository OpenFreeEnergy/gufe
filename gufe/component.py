# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
import hashlib
from typing import Union


class Component(abc.ABC):
    """Base class for members of a ChemicalSystem"""

    # defining this here doesn't actually allow subclasses to inherit, see:
    #  see https://bugs.python.org/issue1549
    def __hash__(self):
        # in python hash(str) isn't deterministic (by design)
        # to provide a stable hash across sessions:
        # generate the string representation of this component
        # (the definition of this string left to individual components)
        s = str(self)

        # pass this through a (fast) stable hash algorithm
        hasher = hashlib.blake2b()
        hasher.update(s.encode('utf-8'))
        dig = hasher.digest()

        # return big endian int representation of the first 8 bytes
        return int.from_bytes(dig[:8], byteorder='big')

    def __lt__(self, other):
        return hash(self) < hash(other)

    @abc.abstractmethod
    def __str__(self):
        ...

    @abc.abstractmethod
    def __eq__(self, other):
        pass

    @abc.abstractmethod
    def to_dict(self) -> dict:
        pass

    @classmethod
    @abc.abstractmethod
    def from_dict(cls, d: dict):
        pass

    @property
    @abc.abstractmethod
    def total_charge(self) -> Union[int, None]:
        """Net formal charge for the Component if defined"""
        ...
