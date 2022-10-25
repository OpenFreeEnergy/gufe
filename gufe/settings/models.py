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
    """
    Settings and modifications we want for all settings classes.
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


class vdWScale(TypedDict):
    """
    See the `SMIRNOFF vdW specification <https://openforcefield.github.io/standards/standards/smirnoff/#vdw>`_
    for more details.
    """
    scale12: float = 0.0
    scale13: float = 0.0
    scale14: float = 0.5
    scale15: float = 1.0


class ElectrostaticScale(TypedDict):
    """
    See the `SMIRNOFF Electrostatics specification <https://openforcefield.github.io/standards/standards/smirnoff/#electrostatics>`_
    for more details.
    """
    scale12: float = 0.0
    scale13: float = 0.0
    scale14: float = 0.833333
    scale15: float = 1.0


class VdWSettings(SettingsBaseModel):
    """
    Settings for van der Waals force.

    See the `SMIRNOFF vdW specification <https://openforcefield.github.io/standards/standards/smirnoff/#vdw>`_
    for more details.
    """

    combining_rules: Literal["Lorentz-Berthelot"]
    potential: Literal["Lennard-Jones-12-6"]
    scale: vdWScale
    long_range_dispersion: Literal["isotropic"]
    cutoff: FloatQuantity["nanometer"] = Field(..., description="Default units are nanometers")
    switch_width: FloatQuantity["nanometer"]
    method: Literal["cutoff"]


class ElectrostaticSettings(SettingsBaseModel):
    """
    Settings for electrostatics
    """

    # Tricky since this allows functions https://openforcefield.github.io/standards/standards/smirnoff/#electrostatics
    periodic_potential: Union[Literal["Ewald3D-ConductingBoundary"], str]
    nonperiodic_potential: Union[Literal["Coulomb"], str]
    exception_potential: Union[Literal["Coulomb"], str]
    scale: ElectrostaticScale
    cutoff: FloatQuantity["nanometer"]  # opt
    switch_width: FloatQuantity["nanometer"]  # opt
    solvent_dielectric: Optional[float]


class GBSASettings(SettingsBaseModel):
    """ "
    Settings for Generalized-Born surface area (GBSA) implicit solvent parameters
    """

    gb_model: str
    solvent_dielectric: float
    solute_dielectric: float


class ForcefieldSettings(SettingsBaseModel):
    """
    Settings for the forcefield
    """

    # Metadata
    date: Optional[date]
    author: Optional[str]

    # These should also allow None
    vdW: VdWSettings
    electrostatics: ElectrostaticSettings
    gbsa: GBSASettings


class ThermoSettings(SettingsBaseModel):
    """
    Settings for thermodynamic parameters
    """

    temperature: FloatQuantity["kelvin"]
    pressure: FloatQuantity["standard_atmosphere"]
    volume: FloatQuantity["nm**3"]
    ph: Union[PositiveFloat, None]
    redox_potential: Optional[float]

    def _validator(self):
        ...


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
    protocol_settings: ProtocolSettings
