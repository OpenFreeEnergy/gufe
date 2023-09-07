import pytest
import copy

from gufe.tokenization import (
    GufeTokenizable,
    new_key_added,
    old_key_removed,
    key_renamed,
    nested_key_moved,
    from_dict,
)

from gufe.tests.test_tokenization import GufeTokenizableTestsMixin


class _DefaultBase(GufeTokenizable):
    """Convenience class to avoid rewriting these methods"""
    @classmethod
    def _from_dict(cls, dct):
        return cls(**dct)

    @classmethod
    def _defaults(cls):
        return super()._defaults()


# this represents an "original" object with  fields `foo` and `bar`
_SERIALIZED_OLD = {
    '__module__': None,  # define in each test
    '__qualname__': None, # define in each test
    'foo': "foo",
    'bar': "bar",
    ':version:': 1,
}


class KeyAdded(_DefaultBase):
    """Add key ``qux`` to the object's dict"""
    def __init__(self, foo, bar, qux=10):
        self.foo = foo
        self.bar = bar
        self.qux = qux

    @classmethod
    def serialization_migration(cls, dct, version):
        if version == 1:
            dct = new_key_added(dct, 'qux', 10)

        return dct

    def _to_dict(self):
        return {"foo": self.foo, "bar": self.bar, "qux": self.qux}


class KeyRemoved(_DefaultBase):
    """Remove key ``bar`` from the object's dict"""
    def __init__(self, foo):
        self.foo = foo

    @classmethod
    def serialization_migration(cls, dct, version):
        if version == 1:
            dct = old_key_removed(dct, "bar", should_warn=True)

        return dct

    def _to_dict(self):
        return {"foo": self.foo}


class KeyRenamed(_DefaultBase):
    """Rename key ``bar`` to ``baz`` in the object's dict"""
    def __init__(self, foo, baz):
        self.foo = foo
        self.baz = baz

    @classmethod
    def serialization_migration(cls, dct, version):
        if version == 1:
            dct = key_renamed(dct, "bar", "baz")

        return dct

    def _to_dict(self):
        return {"foo": self.foo, "baz": self.baz}


class MigrationTester:#(GufeTokenizableTestsMixin):
    cls = None
    """Class to be tested"""
    input_dict = None
    """Initial input dict (except class name info)"""
    kwargs = None
    """kwargs to create an equivalent object from scratch"""

    @property
    def instance(self):
        return self.cls(**self.kwargs)

    def _prep_dct(self, dct):
        dct = copy.deepcopy(self.input_dict)
        dct['__module__'] = self.cls.__module__
        dct['__qualname__'] = self.cls.__qualname__
        return dct

    def test_serialization_migration(self):
        # in these examples, self.kwargs is the same as the output of
        # serialization_migration (not necessarily true for all classes)
        dct = self._prep_dct(self.input_dict)
        del dct['__module__']
        del dct['__qualname__']
        version = dct.pop(':version:')
        assert self.cls.serialization_migration(dct, version) == self.kwargs

    def test_migration(self):
        dct = self._prep_dct(self.input_dict)
        reconstructed = from_dict(dct)
        expected = self.instance
        assert expected == reconstructed

class TestKeyAdded(MigrationTester):
    cls = KeyAdded
    input_dict = _SERIALIZED_OLD
    kwargs = {"foo": "foo", "bar": "bar", "qux": 10}


class TestKeyRemoved(MigrationTester):
    cls = KeyRemoved
    input_dict = _SERIALIZED_OLD
    kwargs = {"foo": "foo"}


class TestKeyRenamed(MigrationTester):
    cls = KeyRenamed
    input_dict = _SERIALIZED_OLD
    kwargs = {"foo": "foo", "baz": "bar"}


# for some reason, we'll move the child from belonging to the son to
# belonging to the daughter (some sort of family issues, idk)
_SERIALIZED_NESTED_OLD = {
    "__module__": ...,
    "__qualname__": ...,
    ":version:": 1,
    "settings": {
        "son": {
            "son_child": 10
        },
        "daughter": {}
    }
}

from pydantic import BaseModel

class SonSettings(BaseModel):
    """v2 model is empty"""

class DaughterSettings(BaseModel):
    """v2 model has child; v1 would not"""
    daughter_child: int

class GrandparentSettings(BaseModel):
    son: SonSettings
    daughter: DaughterSettings

class Grandparent(_DefaultBase):
    def __init__(self, settings: GrandparentSettings):
        self.settings = settings

    def _to_dict(self):
        return {'settings': self.settings.dict()}

    @classmethod
    def _from_dict(cls, dct):
        settings = GrandparentSettings.parse_obj(dct['settings'])
        return cls(settings=settings)

    @classmethod
    @property
    def _version(cls):
        return 2

    @classmethod
    def serialization_migration(cls, dct, version):
        if version == 1:
            dct = nested_key_moved(
                dct,
                old_name="settings.son.son_child",
                new_name="settings.daughter.daughter_child"
            )

        return dct

class TestNestedKeyMoved(MigrationTester):
    cls = Grandparent
    input_dict = _SERIALIZED_NESTED_OLD
    kwargs = {
        'settings': {'son': {}, 'daughter': {'daughter_child': 10}}
    }

    @property
    def instance(self):
        return self.cls(GrandparentSettings(
            son=SonSettings(),
            daughter=DaughterSettings(daughter_child=10)
        ))
