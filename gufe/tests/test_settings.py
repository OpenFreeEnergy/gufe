"""
Tests for convenience functions used to load settings
"""
import pytest

from gufe.settings.settings import get_settings


def test_settings_from_offxml(offxml_settings_path):
    get_settings(offxml_settings_path)


def test_settings_from_json(all_settings_path):
    with pytest.warns(UserWarning):
        get_settings(all_settings_path)
