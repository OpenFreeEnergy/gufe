"""
Pydantic models used for storing settings.
"""

import abc
from datetime import date
from typing import Literal, Optional, Union

from openff.models.models import DefaultModel
from openff.models.types import FloatQuantity
from pydantic import Extra, FilePath, PositiveFloat, Field
from typing_extensions import TypedDict


class SettingsBaseModel(DefaultModel):
    """Settings and modifications we want for all settings classes.
    """

    class Config:
        extra = Extra.forbid
        # Immutability in python is never strict.
        # If developers are determined/stupid
        # they can always modify a so-called "immutable" object.
        # https://pydantic-docs.helpmanual.io/usage/models/#faux-immutability
        allow_mutation = False
        arbitrary_types_allowed = False
        smart_union = True


class vdWScale(SettingsBaseModel):
    """
    See the `SMIRNOFF vdW specification <https://openforcefield.github.io/standards/standards/smirnoff/#vdw>`_
    for more details.
    """
    scale12: float = 0.0
    scale13: float = 0.0
    scale14: float = 0.5
    scale15: float = 1.0


class ElectrostaticScale(SettingsBaseModel):
    """
    See the `SMIRNOFF Electrostatics specification <https://openforcefield.github.io/standards/standards/smirnoff/#electrostatics>`_
    for more details.
    """
    scale12: float = 0.0
    scale13: float = 0.0
    scale14: float = 0.833333
    scale15: float = 1.0


class VdWSettings(SettingsBaseModel):
    """Settings for van der Waals force.

    See the `SMIRNOFF vdW specification <https://openforcefield.github.io/standards/standards/smirnoff/#vdw>`_
    for more details.
    """

    combining_rules: Literal["Lorentz-Berthelot"]
    potential: Literal["Lennard-Jones-12-6"]
    scale: vdWScale
    long_range_dispersion: Literal["isotropic"]
    cutoff: FloatQuantity["angstrom"] = Field(9.0, description="Default units are angstroms")
    switch_width: FloatQuantity["angstrom"] = Field(1.0, description="Default units are angstroms")
    method: Literal["cutoff"]


class ElectrostaticSettings(SettingsBaseModel):
    """Settings for electrostatics.

    See the `SMIRNOFF Electrostatics specification <https://openforcefield.github.io/standards/standards/smirnoff/#electrostatics>`_
    for more details.

    """

    # Tricky since this allows functions https://openforcefield.github.io/standards/standards/smirnoff/#electrostatics
    periodic_potential: Union[Literal["Ewald3D-ConductingBoundary"], str]
    nonperiodic_potential: Union[Literal["Coulomb"], str]
    exception_potential: Union[Literal["Coulomb"], str]
    scale: ElectrostaticScale
    cutoff: FloatQuantity["angstrom"] = Field(None, description="Default units are in angstroms")
    switch_width: FloatQuantity["angstrom"] = Field(None, description="Default units are in angstroms")
    solvent_dielectric: Optional[float]


class GBSASettings(SettingsBaseModel):
    """Settings for Generalized-Born surface area (GBSA) implicit solvent parameters.

    See the `SMIRNOFF GBSA specification <https://openforcefield.github.io/standards/standards/smirnoff/#gbsa>`_
    for more details.
    """

    gb_model: str = "OBC1"
    solvent_dielectric: float = 78.5
    solute_dielectric: float = 1


class ForcefieldSettings(SettingsBaseModel):
    """Settings for the forcefield.
    """

    # Metadata
    date: Optional[date]
    author: Optional[str]

    # These should also allow None
    vdW: VdWSettings
    electrostatics: ElectrostaticSettings
    gbsa: GBSASettings


class ThermoSettings(SettingsBaseModel):
    """Settings for thermodynamic parameters.

    No checking is done to ensure a valid thermodynamic ensemble is possible.
    """

    temperature: FloatQuantity["kelvin"] = Field(None, description="Simulation temperature, default units kelvin")
    pressure: FloatQuantity["standard_atmosphere"] = Field(None, description="Simulation pressure, default units standard atmosphere (atm)")
    volume: FloatQuantity["nm**3"] = Field(None, description="Simulation volume, default units are cubic nanometres")
    ph: Union[PositiveFloat, None] = Field(None, description="Simulation pH, only used if MD engine supports constant pH")
    redox_potential: Optional[float] = Field(None, description="Simulation redox potential, only used if MD engine supports constant E_h")


class ProtocolSettings(SettingsBaseModel, abc.ABC):
    """Protocol-specific settings; this is a base class for protocol
    developers to use for building any settings not included elsewhere.

    """
    ...


class Settings(SettingsBaseModel):
    """
    Container for all settings needed by a protocol
    """

    # symvar? calver?
    settings_version: int
    forcefield_file: Union[FilePath, str]
    forcefield_settings: ForcefieldSettings
    thermo_settings: ThermoSettings
    protocol_settings: Union[ProtocolSettings, None]
