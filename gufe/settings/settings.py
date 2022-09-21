"""
Convenience functions for working with Settings
"""

import json
import warnings
from pathlib import Path

from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.utils.exceptions import SMIRNOFFParseError

from .models import (
    ElectrostaticSettings,
    ForcefieldSettings,
    GBSASettings,
    Settings,
    ThermoSettings,
    VdWSettings,
)


def ff_settings_from_offxml(force_field) -> ForcefieldSettings:
    """
    Creates a ForcefieldSettings object from an open forcefield ForceField object.
    """
    vdW = VdWSettings(
        cutoff=force_field.get_parameter_handler("vdW").cutoff.magnitude,
        combining_rules=force_field.get_parameter_handler("vdW").combining_rules,
        potential=force_field.get_parameter_handler("vdW").potential,
        scale={
            "scale12": force_field.get_parameter_handler("vdW").scale12,
            "scale13": force_field.get_parameter_handler("vdW").scale13,
            "scale14": force_field.get_parameter_handler("vdW").scale14,
            "scale15": force_field.get_parameter_handler("vdW").scale15,
        },
        # see https://github.com/openforcefield/standards/pull/38
        long_range_dispersion="isotropic",
        switch_width=force_field.get_parameter_handler("vdW").switch_width.magnitude,
        method=force_field.get_parameter_handler("vdW").method,
    )
    electrostatics = ElectrostaticSettings(
        periodic_potential=force_field.get_parameter_handler(
            "Electrostatics"
        ).periodic_potential,
        nonperiodic_potential=force_field.get_parameter_handler(
            "Electrostatics"
        ).nonperiodic_potential,
        exception_potential=force_field.get_parameter_handler(
            "Electrostatics"
        ).exception_potential,
        scale={
            "scale12": force_field.get_parameter_handler("Electrostatics").scale12,
            "scale13": force_field.get_parameter_handler("Electrostatics").scale13,
            "scale14": force_field.get_parameter_handler("Electrostatics").scale14,
            "scale15": force_field.get_parameter_handler("Electrostatics").scale15,
        },
        cutoff=force_field.get_parameter_handler("Electrostatics").cutoff.magnitude,
        switch_width=force_field.get_parameter_handler(
            "Electrostatics"
        ).switch_width.magnitude,
        solvent_dielectric=force_field.get_parameter_handler(
            "Electrostatics"
        ).solvent_dielectric,  # Units?
    )
    gbsa = GBSASettings(
        gb_model=force_field.get_parameter_handler("GBSA").gb_model,
        solvent_dielectric=force_field.get_parameter_handler("GBSA").solvent_dielectric,
        solute_dielectric=force_field.get_parameter_handler("GBSA").solute_dielectric,
    )
    forcefield_settings = ForcefieldSettings(
        date=force_field.date,
        author=force_field.author,
        vdW=vdW,
        electrostatics=electrostatics,
        gbsa=gbsa,
    )
    return forcefield_settings


def get_settings(settings_file) -> Settings:
    """
    Read settings file to make Settings object.
    Will try and parse the forcefield file first as offxml to grab settings
    from the forcefield using `ff_settings_from_offxml`.
    """
    settings_file_path = Path(settings_file)
    with open(settings_file_path, "r") as fd:
        settings_from_file = json.load(fd)

    ff_file_path = settings_from_file["forcefield_file"]

    try:
        force_field = ForceField(ff_file_path)
        forcefield_settings = ff_settings_from_offxml(force_field)
    except (SMIRNOFFParseError, OSError):
        # Do we catch OSError as well? Or do we make users have the file exist when they make this object?
        warnings.warn(
            f"not an offxml, all settings will be loaded from {settings_file}",
            UserWarning,
        )
        forcefield_settings = ForcefieldSettings.parse_obj(
            settings_from_file["forcefield_settings"]
        )

    thermo_settings = ThermoSettings.parse_obj(settings_from_file["thermo_settings"])
    settings = Settings(
        settings_version=settings_from_file["settings_version"],
        forcefield_settings=forcefield_settings,
        thermo_settings=thermo_settings,
        forcefield_file=ff_file_path,
    )
    return settings
