# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from gufe.settings import SettingsBaseModel
import pytest


class Foo(SettingsBaseModel):
    a: int
    b: int


class Bar(SettingsBaseModel):
    a: int
    c: int
    d: int


class Baz(SettingsBaseModel):
    thing: Foo
    other: Bar
    d: int


@pytest.fixture
def baz_settings():
    return Baz(
        thing=Foo(a=1, b=2),
        other=Bar(a=10, c=20, d=30),
        d=100,
    )


def test_field_find(baz_settings):
    ret = baz_settings.find_field('a')

    assert set(ret) == {('thing', 'a'), ('other', 'a')}

    ret = baz_settings.find_field('d')

    assert set(ret) == {('d',), ('other', 'd')}

    ret = baz_settings.find_field('e')

    assert ret == []


def test_set_field(baz_settings):
    baz_settings.set_field('b', 99)

    assert baz_settings.thing.b == 99


def test_set_field_not_found(baz_settings):
    with pytest.raises(AttributeError, match="Failed to find"):
        baz_settings.set_field('e', 99)


def test_set_field_ambiguous(baz_settings):
    with pytest.raises(ValueError, match="Ambiguous fieldname"):
        baz_settings.set_field('d', 99)
