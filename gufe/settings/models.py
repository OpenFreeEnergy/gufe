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

    .. note::
       No checking is done to ensure a valid thermodynamic ensemble is
       possible.
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


class BaseForcefieldSettings(SettingsBaseModel, abc.ABC):
    """Base class for ForcefieldSettings objects"""
    ...


class OpenMMSystemGeneratorFFSettings(BaseForcefieldSettings):
    """Parameters to set up the forcefield with OpenMM Forcefields

    .. note::
       Right now we just basically just grab what we need for the
       :class:`openmmforcefields.system_generators.SystemGenerator` 
       signature. See the `OpenMMForceField SystemGenerator documentation`_
       for more details.

    
    .. _`OpenMMForceField SystemGenerator documentation`:
       https://github.com/openmm/openmmforcefields#automating-force-field-management-with-systemgenerator
    """
    constraints: Optional[str] = 'hbonds'
    """Constraints to be applied to system.
       One of 'hbonds', 'allbonds', 'hangles' or None, default 'hbonds'"""

    rigid_water: bool = True
    remove_com: bool = False
    hydrogen_mass: float = 4.0
    """Mass to be repartitioned to hydrogens from neighbouring
       heavy atoms (in amu), default 4.0"""

    forcefields: list[str] = [
        "amber/ff14SB.xml",  # ff14SB protein force field
        "amber/tip3p_standard.xml",  # TIP3P and recommended monovalent ion parameters
        "amber/tip3p_HFE_multivalent.xml",  # for divalent ions
        "amber/phosaa10.xml",  # Handles THE TPO
    ]
    """List of force field paths for all components except :class:`SmallMoleculeComponent`s"""

    small_molecule_forcefield: str = "openff-2.0.0"  # other default ideas 'openff-2.0.0', 'gaff-2.11', 'espaloma-0.2.0'
    """Name of the force field to be used for :class:`SmallMoleculeComponent`s"""

    @pydantic.validator('constraints')
    def constraint_check(cls, v):
        allowed = {'hbonds', 'hangles', 'allbonds'}

        if not (v is None or v.lower() in allowed):
            raise ValueError(f"Bad constraints value, use one of {allowed}")

        return v

class Settings(SettingsBaseModel):
    """
    Container for all settings needed by a protocol

    This represents the minimal surface that all settings objects will have.

    Protocols can subclass this to extend this to cater for their additional settings.
    """
    forcefield_settings: BaseForcefieldSettings
    thermo_settings: ThermoSettings

    @classmethod
    def get_defaults(cls):
        return Settings(
            forcefield_settings=OpenMMSystemGeneratorFFSettings(),
            thermo_settings=ThermoSettings(temperature=300 * unit.kelvin),
        )
