"""
Tests for models used in settings.
Different than testing settings, so tests like model schema generation, model
json round trip, and physical unit testing belongs here.
"""

import json

import pytest
from openff.units import unit

from gufe.settings.models import OpenMMSystemGeneratorFFSettings, Settings, ThermoSettings


def test_settings_schema():
    """Settings schema should be stable"""
    expected_schema = {
        "title": "Settings",
        "description": "Container for all settings needed by a protocol\n\nThis represents the minimal surface that all settings objects will have.\n\nProtocols can subclass this to extend this to cater for their additional settings.",
        "type": "object",
        "properties": {
            "forcefield_settings": {"$ref": "#/definitions/BaseForceFieldSettings"},
            "thermo_settings": {"$ref": "#/definitions/ThermoSettings"},
        },
        "required": ["forcefield_settings", "thermo_settings"],
        "additionalProperties": False,
        "definitions": {
            "BaseForceFieldSettings": {
                "title": "BaseForceFieldSettings",
                "description": "Base class for ForceFieldSettings objects",
                "type": "object",
                "properties": {},
                "additionalProperties": False,
            },
            "ThermoSettings": {
                "title": "ThermoSettings",
                "description": "Settings for thermodynamic parameters.\n\n.. note::\n   No checking is done to ensure a valid thermodynamic ensemble is\n   possible.",
                "type": "object",
                "properties": {
                    "temperature": {
                        "title": "Temperature",
                        "description": "Simulation temperature, default units kelvin",
                        "type": "number",
                    },
                    "pressure": {
                        "title": "Pressure",
                        "description": "Simulation pressure, default units standard atmosphere (atm)",
                        "type": "number",
                    },
                    "ph": {"title": "Ph", "description": "Simulation pH", "exclusiveMinimum": 0, "type": "number"},
                    "redox_potential": {
                        "title": "Redox Potential",
                        "description": "Simulation redox potential",
                        "type": "number",
                    },
                },
                "additionalProperties": False,
            },
        },
    }
    schema = Settings.schema()
    assert schema == expected_schema


def test_default_settings():
    my_settings = Settings.get_defaults()
    my_settings.thermo_settings.temperature = 298 * unit.kelvin
    my_settings.model_dump_json()
    json.dumps(my_settings.model_json_schema(mode='serialization'), indent=2)


class TestSettingsValidation:
    @pytest.mark.parametrize(
        "value,valid,expected",
        [
            ("parsnips", False, None),  # shouldn't be allowed
            ("hbonds", True, "hbonds"),
            ("hangles", True, "hangles"),
            ("allbonds", True, "allbonds"),  # allowed options
            ("HBonds", True, "HBonds"),  # check case insensitivity TODO: cast this to lower?
            (None, True, None),
        ],
    )
    def test_openmmff_constraints(self, value, valid, expected):
        if valid:
            s = OpenMMSystemGeneratorFFSettings(constraints=value)
            assert s.constraints == expected
        else:
            with pytest.raises(ValueError):
                _ = OpenMMSystemGeneratorFFSettings(constraints=value)

    @pytest.mark.parametrize(
        "value,valid,expected",
        [
            (1.0 * unit.nanometer, True, 1.0 * unit.nanometer),
            (0 * unit.nanometer, True, 0 * unit.nanometer),
            (1.0, False, None),  # requires a length unit.
            ("1.1 nm", True, 1.1 * unit.nanometer),
            ("1.1", False, None),
            (-1.0 * unit.nanometer, False, None),
            (1.0 * unit.angstrom, True, 1.0 * unit.angstrom),  # TODO: should we convert this to nm?
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
            with pytest.raises(ValueError):
                _ = OpenMMSystemGeneratorFFSettings(nonbonded_cutoff=value)

    @pytest.mark.parametrize(
        "value,valid,expected",
        [
            ("NoCutoff", True, "NoCutoff"),
            ("NOCUTOFF", False, "NOCUTOFF"),
            ("no cutoff", False, None),
            (1.0, False, None),
        ],
    )
    def test_openmmff_nonbonded_method(self, value, valid, expected):
        if valid:
            s = OpenMMSystemGeneratorFFSettings(nonbonded_method=value)
            assert s.nonbonded_method == expected
        else:
            with pytest.raises(ValueError):
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
            s = ThermoSettings(temperature=value)
            assert s.temperature == expected
        else:
            with pytest.raises(ValueError):
                _ = ThermoSettings(temperature=value)

    @pytest.mark.parametrize(
        "value,valid,expected",
        [
            (1.0 * unit.atm, True, 1.0 * unit.atm),
            (1.0, True, 1.0 * unit.atm),
            ("1 atm", True, 1.0 * unit.atm),
            ("1.0", False, None),
        ],
    )
    def test_thermo_pressure(self, value, valid, expected):
        if valid:
            s = ThermoSettings(pressure=value)
            assert s.pressure == expected
        else:
            with pytest.raises(ValueError):
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
            with pytest.raises(ValueError):
                _ = ThermoSettings(ph=value)

    @pytest.mark.parametrize(
        "value,valid,expected",
        [
            (1.0, True, 1.0),
            (None, True, None),
            ("1", True, 1.0),
        ],
    )
    def test_thermo_redox(self, value, valid, expected):
        if valid:
            s = ThermoSettings(redox_potential=value)
            assert s.redox_potential == expected
        else:
            with pytest.raises(ValueError):
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
        s = Settings.get_defaults()
        s2 = s.frozen_copy()

        assert s == s2

    def test_set_subsection(self):
        # check that attempting to set a subsection of settings still respects
        # frozen state of parent object
        s = Settings.get_defaults().frozen_copy()

        assert s.is_frozen

        ts = ThermoSettings(temperature=301 * unit.kelvin)

        with pytest.raises(AttributeError, match="immutable"):
            s.thermo_settings = ts
