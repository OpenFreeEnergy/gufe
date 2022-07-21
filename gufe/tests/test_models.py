"""
Tests for models used in settings.
Different than testing settings, so tests like model schema generation, model
json round trip, and physical unit testing belongs here.
"""

from gufe.settings.models import Settings
import json

def test_model_schema():
    Settings.schema_json(indent=2)


def test_json_round_trip(all_settings_path, tmp_path):
    with open(all_settings_path) as fd:
        settings = Settings.parse_raw(fd.read())

    assert settings == Settings(**settings.dict())

    d = tmp_path / "test"
    d.mkdir()
    with open(d/"settings.json", "w") as fd:
        fd.write(settings.json())

    with open(d/"settings.json") as fd:
        settings_from_file = json.load(fd)

    assert settings == Settings(**settings_from_file)
