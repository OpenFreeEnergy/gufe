# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from typing import Optional, Iterable

from pydantic import BaseModel, validator
from openff.toolkit.utils.serialization import Serializable

from ..chemicalsystem import ChemicalSystem
from ..mapping import Mapping

from .protocolunit import ProtocolUnit, ProtocolDAG


class Protocol(Serializable, abc.ABC):
    """A protocol that implements an alchemical transformation.

    Takes a `ProtocolSettings` object specific to the protocol on init.
    This configures the protocol for repeated execution on `ChemicalSystem`s.

    Attributes
    ----------
    settings : ProtocolSettings

    """
    ...

    def __init__(
            self,
            settings: "ProtocolSettings" = None
        ):
        """

        """
        self._settings = settings

    def to_dict(self) -> dict:
        ...

    @classmethod
    def from_dict(cls, d: dict):
        ...

    def __hash__(self):
        return hash(
            (
                self._settings
            )
        )

    @classmethod
    @abc.abstractmethod
    def get_default_settings(cls):
        """Get the default settings for this protocol.

        These can be modified and passed back in to the class init.

        """
        ...

    @abc.abstractmethod
    def create(self, 
            initial: ChemicalSystem, 
            final: ChemicalSystem,
            mapping: Mapping = None,
            **kwargs
        ) -> ProtocolDAG:
        """Prepare an `ProtocolDAG` with all information required for
        execution.

        A `ProtocolDAG` is composed of `ProtocolUnit`s, with dependencies on
        

        """
        ...

    def extend(self) -> ProtocolDAG:
        ...
