"""
Pydantic models used for storing settings.
"""

from datetime import date
from typing import Literal, Union

from pydantic import BaseModel, Extra, FilePath, PositiveFloat
from typing_extensions import TypedDict


class SettingsBaseModel(BaseModel):
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
        arbitrary_types_allowed = True
        smart_union = True


class vdWScale(TypedDict):
    scale12: float = 0
    scale13: float = 0
    scale14: float = 0.5
    scale15: float = 1.0


class ElectrostaticScale(TypedDict):
    scale12: float = 0
    scale13: float = 0
    scale14: float = 0.833333
    scale15: float = 1.0


class VdWSettings(SettingsBaseModel):
    """
    Settings for van der Waals force
    """

    combining_rules: Literal["Lorentz-Berthelot"] = "Lorentz-Berthelot"
    potential: Literal["Lennard-Jones-12-6"] = "Lennard-Jones-12-6"
    scale: vdWScale = {"scale12": 0, "scale13": 0, "scale14": 0.5, "scale15": 1.0}
    long_range_dispersion: Literal["isotropic"] = "isotropic"
    cutoff: float = 9.0  # angstrom
    switch_width: float = 1.0  # angstrom
    method: Literal["cutoff"] = "cutoff"


class ElectrostaticSettings(SettingsBaseModel):
    """
    Settings for electrostatics
    """

    # Tricky since this allows functions https://openforcefield.github.io/standards/standards/smirnoff/#electrostatics
    periodic_potential: Union[
        Literal["Ewald3D-ConductingBoundary"], str
    ] = "Ewald3D-ConductingBoundary"
    nonperiodic_potential: Union[Literal["Coulomb"], str] = "Coulomb"
    exception_potential: Union[Literal["Coulomb"], str] = "Coulomb"
    scale: ElectrostaticScale = {
        "scale12": 0,
        "scale13": 0,
        "scale14": 0.833333,
        "scale15": 1.0,
    }
    cutoff: float = None
    switch_width: float = None
    solvent_dielectric: float = None


class BondSettings(SettingsBaseModel):
    """
    Settings for harmonic bonds
    """

    # U(r) = (k/2)*(r-length)^2
    potential: Literal["harmonic"] = "harmonic"
    # Might support other methods?
    fractional_bondorder_method: Literal["AM1-Wiberg"] = "AM1-Wiberg"
    fractional_bondorder_interpolation: Literal["linear"] = "linear"


class AngleSettings(SettingsBaseModel):
    """
    Settings for harmonic angles
    """

    # U(r) = (k/2)*(theta-angle)^2
    potential: Literal["harmonic"] = "harmonic"


class ProperTorsionSettings(SettingsBaseModel):
    """
    Settings for proper torsions
    """

    # U = \sum_{i=1}^N k_i * (1 + cos(periodicity_i * phi - phase_i))
    potential: Literal[
        "k*(1+cos(periodicity*theta-phase))"
    ] = "k*(1+cos(periodicity*theta-phase))"
    default_idivf: Union[Literal["auto"], int] = "auto"
    # Might support other methods?
    fractional_bondorder_method: Literal["AM1-Wiberg"] = "AM1-Wiberg"
    fractional_bondorder_interpolation: Literal["linear"] = "linear"


class ImproperTorsionSettings(SettingsBaseModel):
    """
    Settings for improper torsions
    """

    # U = \sum_{i=1}^N k_i * (1 + cos(periodicity_i * phi - phase_i))
    potential: Literal[
        "k*(1+cos(periodicity*theta-phase))"
    ] = "k*(1+cos(periodicity*theta-phase))"
    # or is it potential="k*(1+cos(periodicity*theta-phase))" ?
    default_idivf: Literal["auto"] = "auto"


class GBSASettings(SettingsBaseModel):
    """ "
    Settings for Generalized-Born surface area (GBSA) implicit solvent parameters
    """

    gb_model: Literal["HCT", "OBC1", "OBC2"] = "OBC1"
    solvent_dielectric: float = 78.5
    solute_dielectric: float = 1


class ForcefieldSettings(SettingsBaseModel):
    """
    Settings for the forcefield
    """

    # Metadata
    date: Union[date, None]
    author: Union[str, None]

    # These should also allow None
    vdW: VdWSettings
    electrostatics: ElectrostaticSettings
    bonds: BondSettings
    angles: AngleSettings
    proper_torsions: ProperTorsionSettings
    improper_torsions: ImproperTorsionSettings
    gbsa: GBSASettings

    # Misc
    aromaticity_model: Literal["OEAroModel_MDL"] = "OEAroModel_MDL"


class ThermoSettings(SettingsBaseModel):
    """
    Settings for thermodynamic parameters
    """

    temperature: PositiveFloat  # Θ
    pressure: PositiveFloat  # M L**−1 T**−2
    ph: PositiveFloat = None
    redox_potential: float = None


class ChemicalComposition(SettingsBaseModel):
    # What does this look like?
    # We take in smiles or file?
    # Or do we take in a SmallMoleculeComponent object?
    pass


class Settings(SettingsBaseModel):
    """
    Container for all settings needed by a protocol
    """

    # symvar? calver?
    settings_version: int = 0
    # Rip this out, only return ff settings from function
    forcefield_file: Union[FilePath, str]
    #####
    forcefield_settings: ForcefieldSettings
    thermo_settings: ThermoSettings
