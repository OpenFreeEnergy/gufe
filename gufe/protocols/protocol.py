# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from typing import Optional

from pydantic import BaseModel, validator
from openff.toolkit.utils.serialization import Serializable

from ..chemicalsystem import ChemicalSystem
from ..mapping import AtomMapping


class Protocol(Serializable, abc.ABC):
    """A protocol that implements an alchemical transformation.

    Takes a `ProtocolSettings` object specific to the protocol on init.
    This configures the protocol for repeated execution on `ChemicalSystem`s.

    Attributes
    ----------
    settings : ProtocolSettings
        THIS MAY BE PROMOTED INTO THE PROTOCOL ITSELF

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
    def is_complete(self) -> bool:
        """Check if the results of this workload already exist"""
        ...

    @abc.abstractmethod
    def prepare(self, 
            initial: ChemicalSystem, 
            final: ChemicalSystem,
            mapping: AtomMapping = None,
            **kwargs
        ) -> Iterable[WorkUnit]:
        """Prepare an iterable of WorkUnits with all information required for
        execution.

        A WorkUnit is the computation required for a Transformation; may map
        to one or more simulations depending on the protocol.

        """
        ...

