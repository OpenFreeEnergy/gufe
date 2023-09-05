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
        expected = self.cls(**self.kwargs)
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
