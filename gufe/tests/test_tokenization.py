import pytest
import abc
from unittest import mock
import json

from gufe.tokenization import (
    GufeTokenizable, GufeKey, tokenize, TOKENIZABLE_REGISTRY,
    import_qualname, get_class, TOKENIZABLE_CLASS_REGISTRY, JSON_HANDLER
)


class Leaf(GufeTokenizable):
    def __init__(self, a, b=2):
        self.a = a
        self.b = b

    def _to_dict(self):
        return {"a": self.a}

    @classmethod
    def _from_dict(cls, dct):
        return cls(**dct)

    def __repr__(self):
        return f"Leaf({self.a})"

    def _defaults(self):
        return super()._defaults()


class Leaf2(GufeTokenizable):
    def __init__(self, a, b=2):
        self.a = a
        self.b = b

    def _to_dict(self):
        return {"a": self.a}

    @classmethod
    def _from_dict(cls, dct):
        return cls(**dct)

    def __repr__(self):
        return f"Leaf({self.a})"

    def _defaults(self):
        return super()._defaults()


class Container(GufeTokenizable):
    def __init__(self, obj, lst, dct):
        self.obj = obj
        self.lst = lst
        self.dct = dct

    def _to_dict(self):
        return {'obj': self.obj, 'lst': self.lst, 'dct': self.dct}

    @classmethod
    def _from_dict(cls, dct):
        return cls(**dct)

    def __repr__(self):
        return f"Container({self.obj}, {self.lst}, {self.dct})"

    def _defaults(self):
        return super()._defaults()


class GufeTokenizableTestsMixin(abc.ABC):

    # set this to the `GufeTokenizable` subclass you are testing
    cls: type[GufeTokenizable]
    key: str

    @pytest.fixture
    def instance(self):
        """Define instance to test with here.

        """
        ...

    def teardown(self):
        TOKENIZABLE_REGISTRY.clear()

    def test_to_dict_roundtrip(self, instance):
        ser = instance.to_dict()
        deser = self.cls.from_dict(ser)
        reser = deser.to_dict()

        assert instance == deser
        assert instance is deser

        # not generally true that the dict forms are equal, e.g. if they
        # include `np.nan`s
        #assert ser == reser

    def test_to_dict_roundtrip_clear_registry(self, instance):
        ser = instance.to_dict()
        TOKENIZABLE_REGISTRY.clear()
        deser = self.cls.from_dict(ser)
        reser = deser.to_dict()

        assert instance == deser
        assert instance is not deser

    def test_to_keyed_dict_roundtrip(self, instance):
        ser = instance.to_keyed_dict()
        deser = self.cls.from_keyed_dict(ser)
        reser = deser.to_keyed_dict()

        assert instance == deser
        assert instance is deser

        # not generally true that the dict forms are equal, e.g. if they
        # include `np.nan`s
        #assert ser == reser

    def test_to_shallow_dict_roundtrip(self, instance):
        ser = instance.to_shallow_dict()
        deser = self.cls.from_shallow_dict(ser)
        reser = deser.to_shallow_dict()

        assert instance == deser
        assert instance is deser

        # not generally true that the dict forms are equal, e.g. if they
        # include `np.nan`s
        #assert ser == reser

    def test_key_stable(self, instance):
        assert self.key == instance.key


class TestGufeTokenizable(GufeTokenizableTestsMixin):

    cls = Container
    key = "Container-262ecded6cd03a619b99d667ded94c9e"

    @pytest.fixture
    def instance(self):
        """Define instance to test with here.

        """
        return self.cont

    def setup(self):
        leaf = Leaf("foo")
        bar = Leaf(leaf)

        self.cont = Container(bar, [leaf, 0], {"leaf": leaf, "a": "b"})

        def leaf_dict(a):
            return {'__module__': __name__, '__qualname__': "Leaf", "a": a}

        self.expected_deep = {
            '__qualname__': "Container",
            '__module__': __name__,
            'obj': leaf_dict(leaf_dict("foo")),
            'lst': [leaf_dict("foo"), 0],
            'dct': {"leaf": leaf_dict("foo"), "a": "b"}
        }

        self.expected_shallow = {
            '__qualname__': "Container",
            '__module__': __name__,
            'obj': bar,
            'lst': [leaf, 0],
            'dct': {'leaf': leaf, 'a': 'b'},
        }

        self.expected_keyed = {
            '__qualname__': "Container",
            '__module__': __name__,
            'obj': {":gufe-key:": bar.key},
            'lst': [{":gufe-key:": leaf.key}, 0],
            'dct': {'leaf': {":gufe-key:": leaf.key}, 'a': 'b'}
        }

    def test_to_dict_deep(self):
        assert self.cont.to_dict() == self.expected_deep

    def test_from_dict_deep(self):
        recreated = Container.from_dict(self.expected_deep)
        assert recreated == self.cont
        assert recreated is self.cont

    def test_to_keyed_dict(self):
        assert self.cont.to_keyed_dict() == self.expected_keyed

    def test_from_keyed_dict(self):
        recreated = self.cls.from_keyed_dict(self.expected_keyed)
        assert recreated == self.cont
        assert recreated is self.cont

    def test_to_shallow_dict(self):
        assert self.cont.to_shallow_dict() == self.expected_shallow

    def test_from_shallow_dict(self):
        recreated = self.cls.from_shallow_dict(self.expected_shallow)
        assert recreated == self.cont
        assert recreated is self.cont

        # here we keep the same objects in memory
        assert recreated.obj.a is recreated.lst[0]
        assert recreated.obj.a is recreated.dct['leaf']

    def test_notequal_different_type(self):
        l1 = Leaf(4)
        l2 = Leaf2(4)

        assert l1 != l2



class Outer:
    class Inner:
        pass


@pytest.mark.parametrize('modname, qualname, expected', [
    (__name__, "Outer", Outer),
    (__name__, "Outer.Inner", Outer.Inner),
    ("gufe.tokenization", 'import_qualname', import_qualname),
])
def test_import_qualname(modname, qualname, expected):
    assert import_qualname(modname, qualname) is expected


def test_import_qualname_not_yet_imported():
    # this is specifically to test that something we don't have imported in
    # this module will import correctly
    msg_cls = import_qualname(modname="email.message",
                              qualname="EmailMessage")
    from email.message import EmailMessage
    assert msg_cls is EmailMessage


def test_import_qualname_remappings():
    remappings = {("foo", "Bar.Baz"): (__name__, "Outer.Inner")}
    assert import_qualname("foo", "Bar.Baz", remappings) is Outer.Inner


@pytest.mark.parametrize('modname, qualname', [
    (None, "Outer.Inner"),
    (__name__, None),
])
def test_import_qualname_error_none(modname, qualname):
    with pytest.raises(ValueError, match="cannot be None"):
        import_qualname(modname, qualname)



@pytest.mark.parametrize('cls_reg', [
    {},
    {(__name__, "Outer.Inner"): Outer.Inner},
])
def test_get_class(cls_reg):
    with mock.patch.dict("gufe.tokenization.TOKENIZABLE_CLASS_REGISTRY", cls_reg):
        assert get_class(__name__, "Outer.Inner") is Outer.Inner

def test_path_to_json():
    import pathlib
    p = pathlib.Path("foo/bar")
    ser = json.dumps(p, cls=JSON_HANDLER.encoder)
    deser = json.loads(ser, cls=JSON_HANDLER.decoder)
    reser = json.dumps(deser, cls=JSON_HANDLER.encoder)
    assert ser == reser
    assert p == deser

