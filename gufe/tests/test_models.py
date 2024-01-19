"""
Tests for models used in settings.
Different than testing settings, so tests like model schema generation, model
json round trip, and physical unit testing belongs here.
"""

import json

from openff.units import unit
import pytest

from gufe.settings.models import (
    OpenMMSystemGeneratorFFSettings,
    Settings,
    ThermoSettings,
)


def test_model_schema():
    Settings.schema_json(indent=2)


@pytest.mark.xfail  # issue #125
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


@pytest.mark.parametrize('value,good', [
    ('parsnips', False),  # shouldn't be allowed
    ('hbonds', True), ('hangles', True), ('allbonds', True),  # allowed options
    ('HBonds', True),  # check case insensitivity
])
def test_invalid_constraint(value, good):
    if good:
        s = OpenMMSystemGeneratorFFSettings(constraints=value)
        assert s
    else:
        with pytest.raises(ValueError):
            _ = OpenMMSystemGeneratorFFSettings(constraints=value)


class TestFreezing:
    def test_default_not_frozen(self):
        s = Settings.get_defaults()
        # make a frozen copy to check this doesn't alter the original
        s2 = s.frozen_copy()

        s.thermo_settings.temperature = 199 * unit.kelvin
        assert s.thermo_settings.temperature == 199 * unit.kelvin

    def test_freezing(self):
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
