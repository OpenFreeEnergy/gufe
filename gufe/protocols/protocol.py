# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from typing import Optional

from pydantic import BaseModel, validator
from openff.toolkit.utils.serialization import Serializable

from ..chemicalsystem import ChemicalSystem
from ..mapping import AtomMapping
from .protocol_settings import ProtocolSettings


class Protocol(abc.ABC, Serializable):
    """A protocol that implements an alchemical transformation.

    Takes a `ProtocolSettings` object specific to the protocol on init.
    This configures the protocol for repeated execution on `ChemicalSystem`s.

    The `run` method takes the `start` and `end` `ChemicalSystem`,
    the `AtomMapping`, and any additional keywords specific to the `Protocol`.

    Attributes
    ----------
    settings : 

    """
    ...

    def __init__(
            self,
            *,
            settings: ProtocolSettings = None
        ):
        """

        """
        self._settings = settings

    def to_dict(self) -> dict:
        ...

    @classmethod
    def from_dict(cls, d: dict):
        ...

    @classmethod
    @abc.abstractmethod
    def get_default_settings(cls):
        """Get the default settings for this protocol.

        These can be modified and passed back in to the class init.

        """
        ...

    @abc.abstractmethod
    def is_complete(self) -> bool:
        """Check if the results of this workload already exist"""
        ...

    @abc.abstractmethod
    def run(self, 
            start: ChemicalSystem, 
            end: ChemicalSystem,
            mapping: AtomMapping = None,
            **kwargs
        ) -> bool:
        """Perform this method, returning success.


        """
        ...


class LigandSolventAtomMappedProtocol(Protocol):
    ...


class LigandComplexAtomMappedProtocol(Protocol):
    ...



# want a variant of the above based on engine, probably
# e.g. OpenMM, Gromacs
# since settings will be very different
# not exactly sure how to avoid long names here
