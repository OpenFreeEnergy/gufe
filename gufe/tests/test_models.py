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

@pytest.mark.xfail  # issue #125  TODO: update this now that files are vendored
def test_json_round_trip(all_settings_path, tmp_path):
    with open(all_settings_path) as fd:
        settings = Settings.parse_raw(fd.read())

    assert settings == Settings(**settings.dict())

    d = tmp_path / "test"
    d.mkdir()
    with open(d / "settings.json", "w") as fd:
        fd.write(settings.json())

    with open(d / "settings.json") as fd:
        settings_from_file = json.load(fd)

    assert settings == Settings.parse_raw(settings_from_file)

def test_default_settings():
    my_settings = Settings.get_defaults()
    my_settings.thermo_settings.temperature = 298 * unit.kelvin
    my_settings.json()
    my_settings.schema_json(indent=2)


class TestOpenMMSystemGeneratorFFSettings:
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
    def test_constraints_validation(self, value, valid, expected):
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
            (1.0, True, 1.0 * unit.nanometer),  # should cast float to nanometer
            ("1.1 nm", True, 1.1 * unit.nanometer),
            ("1.1 ", False, None),
            # (1.0 * unit.angstrom, True, 0.100 * unit.nanometer),  # TODO: why does this not work?
            (300 * unit.kelvin, False, None),
            (True, False, None),
            # ("one", False, None),  # TODO: more elegant error handling for this
        ],
    )
    def test_nonbonded_cutoff_validation(self, value, valid, expected):
        if valid:
            s = OpenMMSystemGeneratorFFSettings(nonbonded_cutoff=value)
            assert s.nonbonded_cutoff == expected
        else:
            with pytest.raises(ValueError):
                _ = OpenMMSystemGeneratorFFSettings(nonbonded_cutoff=value)


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
