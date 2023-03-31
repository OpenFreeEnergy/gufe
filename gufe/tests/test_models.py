"""
Tests for models used in settings.
Different than testing settings, so tests like model schema generation, model
json round trip, and physical unit testing belongs here.
"""

import json

from openff.units import unit
import pytest

from gufe.settings.models import Settings, OpenMMSystemGeneratorFFSettings


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
