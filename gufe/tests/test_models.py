"""
Tests for models used in settings.
Different than testing settings, so tests like model schema generation, model
json round trip, and physical unit testing belongs here.
"""

import json

import numpy as np
import pytest
from openff.units import unit
from openmm import unit as openmm_unit

from gufe.settings import SettingsBaseModel
from gufe.settings.models import OpenMMSystemGeneratorFFSettings, Settings, ThermoSettings
from gufe.settings.types import BoxQuantity


def test_settings_schema():
    """Settings schema should be stable"""
    expected_schema = {
        "$defs": {
            "BaseForceFieldSettings": {
                "additionalProperties": False,
                "description": "Base class for ForceFieldSettings objects",
                "properties": {},
                "title": "BaseForceFieldSettings",
                "type": "object",
            },
            "ThermoSettings": {
                "additionalProperties": False,
                "description": "Settings for thermodynamic parameters.\n\n.. note::\n   No checking is done to ensure a valid thermodynamic ensemble is\n   possible.",
                "properties": {
                    "temperature": {
                        "anyOf": [{"type": "number"}, {"type": "null"}],
                        "default": None,
                        "description": "Simulation temperature in kelvin)",
                        "title": "Temperature",
                    },
                    "pressure": {
                        "anyOf": [{"type": "number"}, {"type": "null"}],
                        "default": None,
                        "description": "Simulation pressure in standard atmosphere (atm)",
                        "title": "Pressure",
                    },
                    "ph": {
                        "anyOf": [{"exclusiveMinimum": 0, "type": "number"}, {"type": "null"}],
                        "default": None,
                        "description": "Simulation pH",
                        "title": "Ph",
                    },
                    "redox_potential": {
                        "anyOf": [{"type": "number"}, {"type": "null"}],
                        "default": None,
                        "description": "Simulation redox potential in millivolts (mV).",
                        "title": "Redox Potential",
                    },
                },
                "title": "ThermoSettings",
                "type": "object",
            },
        },
        "additionalProperties": False,
        "description": "Container for all settings needed by a protocol\n\nThis represents the minimal surface that all settings objects will have.\n\nProtocols can subclass this to extend this to cater for their additional settings.",
        "properties": {
            "forcefield_settings": {"$ref": "#/$defs/BaseForceFieldSettings", "title": "Forcefield Settings"},
            "thermo_settings": {"$ref": "#/$defs/ThermoSettings", "title": "Thermo Settings"},
        },
        "required": ["forcefield_settings", "thermo_settings"],
        "title": "Settings",
        "type": "object",
    }
    ser_schema = Settings.model_json_schema(mode="serialization")
    val_schema = Settings.model_json_schema(mode="validation")

    assert ser_schema == expected_schema
    assert val_schema == expected_schema


def test_openmmffsettings_schema():
    expected_schema = {
        "additionalProperties": False,
        "description": "Parameters to set up the force field with OpenMM ForceFields\n\n.. note::\n   Currently, this stores what is needed for the\n   :class:`openmmforcefields.system_generators.SystemGenerator` signature.\n   See the `OpenMMForceField SystemGenerator documentation`_ for more details.\n\n\n.. _`OpenMMForceField SystemGenerator documentation`:\n   https://github.com/openmm/openmmforcefields#automating-force-field-management-with-systemgenerator",
        "properties": {
            "constraints": {
                "anyOf": [{"enum": ["hbonds", "allbonds", "hangles"], "type": "string"}, {"type": "null"}],
                "default": "hbonds",
                "title": "Constraints",
            },
            "rigid_water": {"default": True, "title": "Rigid Water", "type": "boolean"},
            "hydrogen_mass": {"default": 3.0, "title": "Hydrogen Mass", "type": "number"},
            "forcefields": {
                "default": [
                    "amber/ff14SB.xml",
                    "amber/tip3p_standard.xml",
                    "amber/tip3p_HFE_multivalent.xml",
                    "amber/phosaa10.xml",
                ],
                "items": {"type": "string"},
                "title": "Forcefields",
                "type": "array",
            },
            "small_molecule_forcefield": {
                "default": "openff-2.2.1",
                "title": "Small Molecule Forcefield",
                "type": "string",
            },
            "nonbonded_method": {"default": "PME", "title": "Nonbonded Method", "type": "string"},
            "nonbonded_cutoff": {
                "description": "Cutoff value for short range nonbonded interactions.",
                "title": "Nonbonded Cutoff",
                "type": "number",
            },
        },
        "title": "OpenMMSystemGeneratorFFSettings",
        "type": "object",
    }
    ser_schema = OpenMMSystemGeneratorFFSettings.model_json_schema(mode="serialization")
    val_schema = OpenMMSystemGeneratorFFSettings.model_json_schema(mode="validation")
    assert ser_schema == expected_schema
    assert val_schema == expected_schema


def test_default_settings():
    my_settings = Settings.get_defaults()
    my_settings.thermo_settings.temperature = 298 * unit.kelvin
    my_settings.model_dump_json()
    json.dumps(my_settings.model_json_schema(mode="serialization"), indent=2)


class TestSettingsValidation:
    @pytest.mark.parametrize(
        "value,valid,expected",
        [
            ("Parsnips", False, None),  # shouldn't be allowed
            (1.0, False, None),  # shouldn't be allowed
            ("hbonds", True, "hbonds"),
            ("hangles", True, "hangles"),
            ("allbonds", True, "allbonds"),  # allowed options
            ("HBonds", True, "hbonds"),  # check case insensitivity
            (None, True, None),
        ],
    )
    def test_openmmff_constraints(self, value, valid, expected):
        if valid:
            s = OpenMMSystemGeneratorFFSettings(constraints=value)
            assert s.constraints == expected
        else:
            with pytest.raises(ValueError, match="Input should be 'hbonds', 'allbonds' or 'hangles'"):
                _ = OpenMMSystemGeneratorFFSettings(constraints=value)

    @pytest.mark.parametrize(
        "value,valid,expected",
        [
            (1.0 * unit.nanometer, True, 1.0 * unit.nanometer),
            (openmm_unit.Quantity(2.0, openmm_unit.nanometer), True, 2.0 * unit.nanometer),
            ({"val": 1.0, "unit": unit.nanometer}, True, 1.0 * unit.nanometer),
            (1.0, False, None),  # requires a length unit.
            ("1.1 nm", True, 1.1 * unit.nanometer),
            ("1.1", False, None),
            # NOTE: this is not precisely equal for smaller values due to pint unit floating point precision
            (100.0 * unit.angstrom, True, 10.0 * unit.nanometer),
            (300 * unit.kelvin, False, None),
            (True, False, None),
            (None, False, None),
            # ("one", False, None),  # TODO: more elegant error handling for this
        ],
    )
    def test_openmmff_nonbonded_cutoff(self, value, valid, expected):
        if valid:
            s = OpenMMSystemGeneratorFFSettings(nonbonded_cutoff=value)
            assert s.nonbonded_cutoff == expected
        else:
            with pytest.raises(ValueError, match="nonbonded_cutoff"):
                _ = OpenMMSystemGeneratorFFSettings(nonbonded_cutoff=value)

    @pytest.mark.parametrize(
        "value,valid,expected",
        [
            (0 * unit.nanometer, True, 0 * unit.nanometer),
            (-1.0 * unit.nanometer, False, None),
        ],
    )
    def test_negative_cutoff_error(self, value, valid, expected):
        if valid:
            s = OpenMMSystemGeneratorFFSettings(nonbonded_cutoff=value)
            assert s.nonbonded_cutoff == expected
        else:
            with pytest.raises(ValueError, match=" Input should be greater than or equal to 0"):
                _ = OpenMMSystemGeneratorFFSettings(nonbonded_cutoff=value)

    @pytest.mark.parametrize(
        "value,valid,expected",
        [
            ("NoCutoff", True, "NoCutoff"),
            ("NOCUTOFF", True, "NOCUTOFF"),
            ("no cutoff", False, None),
            (1.0, False, None),
        ],
    )
    def test_openmmff_nonbonded_method(self, value, valid, expected):
        if valid:
            s = OpenMMSystemGeneratorFFSettings(nonbonded_method=value)
            assert s.nonbonded_method == expected
        else:
            with pytest.raises(ValueError, match="nonbonded_method"):
                _ = OpenMMSystemGeneratorFFSettings(nonbonded_method=value)

    @pytest.mark.parametrize(
        "value,valid,expected",
        [
            (298 * unit.kelvin, True, 298 * unit.kelvin),
            ("298 kelvin", True, 298 * unit.kelvin),
            (298, False, None),  # requires units
            ("298", False, None),
            (298 * unit.angstrom, False, None),
        ],
    )
    def test_thermo_temperature(self, value, valid, expected):
        if valid:
            settings = ThermoSettings(temperature=value)
            assert settings.temperature == expected
        else:
            with pytest.raises(ValueError, match="temperature"):
                _ = ThermoSettings(temperature=value)

    @pytest.mark.parametrize(
        "value,valid,expected",
        [
            (1.0 * unit.atm, True, 1.0 * unit.atm),
            (1.0, False, None),  # require units
            ("1 atm", True, 1.0 * unit.atm),
            ("1.0", False, None),
        ],
    )
    def test_thermo_pressure(self, value, valid, expected):
        if valid:
            s = ThermoSettings(pressure=value)
            assert s.pressure == expected
        else:
            with pytest.raises(ValueError, match="pressure"):
                _ = ThermoSettings(pressure=value)

    @pytest.mark.parametrize(
        "value,valid,expected",
        [
            (1.0, True, 1.0),
            (1, True, 1.0),
            ("1 ph", False, None),
            (1.0 * unit.atm, False, None),
        ],
    )
    def test_thermo_ph(self, value, valid, expected):
        if valid:
            s = ThermoSettings(ph=value)
            assert s.ph == expected
        else:
            with pytest.raises(ValueError, match="ph"):
                _ = ThermoSettings(ph=value)

    @pytest.mark.parametrize(
        "value,valid,expected",
        [
            (None, True, None),
            (1 * unit.mV, True, 1 * unit.mV),
            ("1.0 mV", True, 1 * unit.mV),
            ("0.001 volts", True, 1 * unit.mV),
            (0.001 * unit.volt, True, 1 * unit.mV),
            (0.001 * unit.nanometer, False, None),
            ("0.001 nm", False, None),
        ],
    )
    def test_thermo_redox(self, value, valid, expected):
        if valid:
            s = ThermoSettings(redox_potential=value)
            assert s.redox_potential == expected
        else:
            with pytest.raises(ValueError, match="redox"):
                _ = ThermoSettings(redox_potential=value)


class TestFreezing:
    def test_default_not_frozen(self):
        s = Settings.get_defaults()
        # make a frozen copy to check this doesn't alter the original
        s2 = s.frozen_copy()

        s.thermo_settings.temperature = 199 * unit.kelvin
        assert s.thermo_settings.temperature == 199 * unit.kelvin

    def test_default_frozen(self):
        s = Settings.get_defaults()

        s2 = s.frozen_copy()

        with pytest.raises(AttributeError, match="immutable"):
            s2.thermo_settings.temperature = 199 * unit.kelvin

    def test_unfreezing(self):
        s = Settings.get_defaults()

        s2 = s.frozen_copy()

        with pytest.raises(AttributeError, match="immutable"):
            s2.thermo_settings.temperature = 199 * unit.kelvin

        assert s2.is_frozen

        s3 = s2.unfrozen_copy()

        s3.thermo_settings.temperature = 199 * unit.kelvin
        assert s3.thermo_settings.temperature == 199 * unit.kelvin

    def test_frozen_equality(self):
        # the frozen-ness of Settings doesn't alter its contents
        # therefore a frozen/unfrozen Settings which are otherwise identical
        # should be considered equal
        s1 = Settings.get_defaults()
        s2 = s1.frozen_copy()

        assert s1 == s2

    def test_frozen_equality_changed(self):
        # the frozen-ness of Settings doesn't alter its contents
        # therefore a frozen/unfrozen Settings which are otherwise identical
        # should be considered equal
        s1 = Settings.get_defaults()
        s2 = s1.frozen_copy()
        s1.forcefield_settings.constraints = "allbonds"
        assert s1 != s2

    def test_settings_equality_not_settings(self):
        """check that our custom __eq__ implementation handles non-settings objects"""
        s1 = Settings.get_defaults()
        assert s1 != "not a settings object"

    def test_set_subsection(self):
        # check that attempting to set a subsection of settings still respects
        # frozen state of parent object
        s = Settings.get_defaults().frozen_copy()

        assert s.is_frozen

        ts = ThermoSettings(temperature=301 * unit.kelvin)

        with pytest.raises(AttributeError, match="immutable"):
            s.thermo_settings = ts


class BoxSettingsModel(SettingsBaseModel):
    box_vectors: BoxQuantity


def test_box_quantity_schema():
    expected_schema = {
        "additionalProperties": False,
        "properties": {"box_vectors": {"title": "Box Vectors", "type": "number"}},
        "required": ["box_vectors"],
        "title": "BoxSettingsModel",
        "type": "object",
    }
    ser_schema = BoxSettingsModel.model_json_schema(mode="serialization")
    val_schema = BoxSettingsModel.model_json_schema(mode="validation")
    assert ser_schema == expected_schema
    assert val_schema == expected_schema


@pytest.mark.parametrize(
    "value",
    [
        np.asarray([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]),
        np.asarray([1.0, 1.0, 1.0]),
        [1.0, 1.0, 1.0],
        [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]] * unit.angstrom,
    ],
)
def test_valid_box_quantity(value):
    box_settings = BoxSettingsModel(box_vectors=value)
    assert box_settings.box_vectors.units == unit.nanometer


# TODO: improve this error handling
# @pytest.mark.parametrize(
#     "value,err_msg",
#     [
#         ("a string", None),
#         (1.0*unit.nanometer, None),
#         (1.0, None),
#     ],
# )
# def test_invalid_box_quantity(value, err_msg):
#     with pytest.raises(RuntimeError):
#         BoxSettingsModel(box_vectors=value)
