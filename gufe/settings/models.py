# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""
Pydantic models used for storing settings.
"""

import abc
from typing import Optional, Union

import pydantic
from openff.models.models import DefaultModel
from openff.models.types import FloatQuantity
from openff.units import unit
from pydantic import Extra, Field, PositiveFloat


class SettingsBaseModel(DefaultModel):
    """Settings and modifications we want for all settings classes."""

    class Config:
        extra = Extra.forbid
        arbitrary_types_allowed = False
        smart_union = True


class ThermoSettings(SettingsBaseModel):
    """Settings for thermodynamic parameters.

    No checking is done to ensure a valid thermodynamic ensemble is possible.
    """

    temperature: FloatQuantity["kelvin"] = Field(
        None, description="Simulation temperature, default units kelvin"
    )
    pressure: FloatQuantity["standard_atmosphere"] = Field(
        None, description="Simulation pressure, default units standard atmosphere (atm)"
    )
    ph: Union[PositiveFloat, None] = Field(None, description="Simulation pH")
    redox_potential: Optional[float] = Field(
        None, description="Simulation redox potential"
    )


class ProtocolSettings(SettingsBaseModel, abc.ABC):
    """Protocol-specific settings; this is a base class for protocol
    developers to use for building any settings not included elsewhere.

    """

    ...


class BaseForcefieldSettings(SettingsBaseModel, abc.ABC):
    """Base class for ForcefieldSettings objects"""

    ...


class OpenMMSystemGeneratorFFSettings(BaseForcefieldSettings):
    """Parameters to set up the forcefield with OpenMM Forcefields

    Attributes
    ----------
    constraints : Optional[str]
      one of 'hbonds', 'allbonds', 'hangles' or None, default 'hbonds'
    rigid_water : bool
      default True
    remove_com : bool
       default False
    hydrogen_mass : float
       the repartitioned hydrogen mass (in amu), default 4.0
    forcefields : list[str]

    small_molecule_forcefield : str

    Right now we just basically just grab what we need for the systemgenerator signature
    https://github.com/openmm/openmmforcefields#automating-force-field-management-with-systemgenerator
    """
    constraints: Optional[str] = 'hbonds'
    rigid_water: bool = True
    remove_com: bool = False
    hydrogen_mass: float

    forcefields: list[str] = [
        "amber/ff14SB.xml",  # ff14SB protein force field
        "amber/tip3p_standard.xml",  # TIP3P and recommended monovalent ion parameters
        "amber/tip3p_HFE_multivalent.xml",  # for divalent ions
        "amber/phosaa10.xml",  # Handles THE TPO
    ]
    small_molecule_forcefield: str = "openff-2.0.0"  # other default ideas 'openff-2.0.0', 'gaff-2.11', 'espaloma-0.2.0'

    @pydantic.validator
    def constraint_check(cls, v):
        allowed = {'hbonds', 'hangles', 'allbonds'}

        if not (v is None or v.lower() in allowed):
            raise ValueError(f"Bad constraints value, use one of {allowed}")

        return v

class Settings(SettingsBaseModel):
    """
    Container for all settings needed by a protocol
    """

    # symvar? calver?
    settings_version: int = 0
    forcefield_settings: OpenMMSystemGeneratorFFSettings
    thermo_settings: ThermoSettings
    protocol_settings: Union[ProtocolSettings, None]

    @classmethod
    def get_defaults(cls):
        return Settings(
            forcefield_settings=OpenMMSystemGeneratorFFSettings(),
            thermo_settings=ThermoSettings(temperature=300 * unit.kelvin),
        )
