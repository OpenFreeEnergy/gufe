import pytest
import abc
import datetime
import logging
import io
from unittest import mock
import json
from typing import Optional

from gufe.tokenization import (
    GufeTokenizable, GufeKey, tokenize, TOKENIZABLE_REGISTRY,
    import_qualname, get_class, TOKENIZABLE_CLASS_REGISTRY, JSON_HANDLER,
    get_all_gufe_objs, gufe_to_digraph, gufe_objects_from_shallow_dict,
    KeyedChain,
)


class Leaf(GufeTokenizable):
    def __init__(self, a, b=2):
        self.logger.info("no key defined!")
        self.a = a
        self.b = b
        self.logger.info(f"{a=}")
        self.logger.debug(f"{b=}")

    def _to_dict(self):
        return {"a": self.a, "b": self.b}

    @classmethod
    def _from_dict(cls, dct):
        return cls(**dct)

    def __repr__(self):
        return f"Leaf({self.a}, {self.b})"

    @classmethod
    def _defaults(cls):
        return super()._defaults()


class Leaf2(GufeTokenizable):
    def __init__(self, a, b=2):
        self.a = a
        self.b = b

    def _to_dict(self):
        return {"a": self.a, "b": self.b}

    @classmethod
    def _from_dict(cls, dct):
        return cls(**dct)

    def __repr__(self):
        return f"Leaf({self.a}, {self.b})"

    @classmethod
    def _defaults(cls):
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

    @classmethod
    def _defaults(cls):
        return super()._defaults()


class GufeTokenizableTestsMixin(abc.ABC):

    # set this to the `GufeTokenizable` subclass you are testing
    cls: type[GufeTokenizable]
    key: Optional[str]
    repr: Optional[str]

    @pytest.fixture
    def instance(self):
        """Define instance to test with here.

        """
        ...

    def test_to_dict_roundtrip(self, instance):
        ser = instance.to_dict()
        deser = self.cls.from_dict(ser)
        reser = deser.to_dict()

        assert instance == deser
        assert instance is deser

        # not generally true that the dict forms are equal, e.g. if they
        # include `np.nan`s
        #assert ser == reser

    @pytest.mark.skip
    def test_to_dict_roundtrip_clear_registry(self, instance):
        ser = instance.to_dict()
        patch_loc = "gufe.tokenization.TOKENIZABLE_REGISTRY"
        with mock.patch.dict(patch_loc, {}, clear=True):
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

    def test_repr(self, instance):
        if self.repr is None:
            # nondeterministic reprs
            assert isinstance(repr(instance), str)
        else:
            assert repr(instance) == self.repr


class TestGufeTokenizable(GufeTokenizableTestsMixin):

    cls = Container
    key = "Container-3fcec08974fbbd0371fed8a185628b70"
    repr = "Container(Leaf(Leaf(foo, 2), 2), [Leaf(foo, 2), 0], {'leaf': Leaf(foo, 2), 'a': 'b'})"

    @pytest.fixture
    def instance(self):
        """Define instance to test with here.

        """
        return self.cont

    def setup_method(self):
        leaf = Leaf("foo")
        bar = Leaf(leaf)

        self.cont = Container(bar, [leaf, 0], {"leaf": leaf, "a": "b"})

        def leaf_dict(a):
            return {'__module__': __name__, '__qualname__': "Leaf", "a": a,
                    "b": 2, ':version:': 1}

        self.expected_deep = {
            '__qualname__': "Container",
            '__module__': __name__,
            'obj': leaf_dict(leaf_dict("foo")),
            'lst': [leaf_dict("foo"), 0],
            'dct': {"leaf": leaf_dict("foo"), "a": "b"},
            ':version:': 1,
        }

        self.expected_shallow = {
            '__qualname__': "Container",
            '__module__': __name__,
            'obj': bar,
            'lst': [leaf, 0],
            'dct': {'leaf': leaf, 'a': 'b'},
            ':version:': 1,
        }

        self.expected_keyed = {
            '__qualname__': "Container",
            '__module__': __name__,
            'obj': {":gufe-key:": bar.key},
            'lst': [{":gufe-key:": leaf.key}, 0],
            'dct': {'leaf': {":gufe-key:": leaf.key}, 'a': 'b'},
            ':version:': 1,
        }

    def test_set_key(self):
        leaf = Leaf("test-set-key")
        key = leaf.key
        patch_loc = "gufe.tokenization.TOKENIZABLE_REGISTRY"
        registry = dict(TOKENIZABLE_REGISTRY)
        with mock.patch.dict(patch_loc, registry, clear=True):
            # import pdb; pdb.set_trace()
            leaf._set_key("qux")
            assert leaf.key == "qux"
            assert TOKENIZABLE_REGISTRY["qux"] is leaf
            assert key not in TOKENIZABLE_REGISTRY

        assert TOKENIZABLE_REGISTRY[key] is leaf
        assert "qux" not in TOKENIZABLE_REGISTRY

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

    def test_copy_with_replacements(self):
        l1 = Leaf(4)
        l2 = l1.copy_with_replacements(b=4)
        assert l1 != l2
        assert l1.a == l2.a
        assert l1.b != l2.b

    def test_copy_with_replacements_no_arguments(self):
        # with no arguments, copy_with_replacements returns the same object
        l1 = Leaf(4)
        l2 = l1.copy_with_replacements()
        assert l1 == l2
        assert l1 is l2

    def test_copy_with_replacements_invalid(self):
        l1 = Leaf(4)
        with pytest.raises(TypeError, match="Invalid"):
            _ = l1.copy_with_replacements(foo=10)

    @pytest.mark.parametrize('level', ["DEBUG", "INFO", "CRITICAL"])
    def test_logging(self, level):
        stream = io.StringIO()
        handler = logging.StreamHandler(stream)
        fmt = logging.Formatter(
            "%(name)s - %(gufekey)s - %(levelname)s - %(message)s"
        )
        name = "gufekey.gufe.tests.test_tokenization.Leaf"
        logger = logging.getLogger(name)
        logger.setLevel(getattr(logging, level))
        handler.setFormatter(fmt)
        logger.addHandler(handler)

        leaf = Leaf(10)

        results = stream.getvalue()
        key = leaf.key.split('-')[-1]

        initial_log = f"{name} - UNKNOWN - INFO - no key defined!\n"
        info_log = f"{name} - {key} - INFO - a=10\n"
        debug_log = f"{name} - {key} - DEBUG - b=2\n"

        expected = ""
        if level in {"DEBUG", "INFO"}:
            expected += initial_log + info_log
        if level == "DEBUG":
            expected += debug_log

        assert results == expected


def test_get_all_gufe_objs():
    leaf = Leaf("foo")
    bar = Leaf(leaf)
    cont = Container(bar, [leaf, 0], {"leaf": leaf, "a": "b"})
    all_objs = get_all_gufe_objs(cont)
    assert all_objs == {cont, bar, leaf}


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


class TestGufeKey:
    def test_to_dict(self):
        k = GufeKey('foo-bar')

        assert k.to_dict() == {':gufe-key:': 'foo-bar'}

    def test_prefix(self):
        k = GufeKey('foo-bar')

        assert k.prefix == 'foo'

    def test_token(self):
        k = GufeKey('foo-bar')

        assert k.token == 'bar'


def test_gufe_to_digraph(solvated_complex):
    graph = gufe_to_digraph(solvated_complex)

    connected_objects = gufe_objects_from_shallow_dict(
        solvated_complex.to_shallow_dict()
    )

    assert len(graph.nodes) == 4
    assert len(graph.edges) == 3

    for node_a, node_b in graph.edges:
        assert node_b in connected_objects
        assert node_a is solvated_complex


def test_gufe_objects_from_shallow_dict(solvated_complex):
    shallow_dict = solvated_complex.to_shallow_dict()
    gufe_objects = set(gufe_objects_from_shallow_dict(shallow_dict))

    assert len(gufe_objects) == 3

    for gufe_object in gufe_objects:
        assert gufe_object in solvated_complex.components.values()


class TestKeyedChain:

    def test_from_gufe(self, benzene_variants_star_map):
        contained_objects = list(get_all_gufe_objs(benzene_variants_star_map))
        expected_len = len(contained_objects)

        kc = KeyedChain.from_gufe(benzene_variants_star_map)

        assert len(kc) == expected_len

        original_keys = [obj.key for obj in contained_objects]
        original_keyed_dicts = [
            obj.to_keyed_dict() for obj in contained_objects
        ]

        kc_gufe_keys = set(kc.gufe_keys())
        kc_shallow_dicts = list(kc.keyed_dicts())

        for key, keyed_dict in zip(original_keys, original_keyed_dicts):
            assert key in kc_gufe_keys
            assert keyed_dict in kc_shallow_dicts

    def test_to_gufe(self, benzene_variants_star_map):
        kc = KeyedChain.from_gufe(benzene_variants_star_map)
        assert hash(kc.to_gufe()) == hash(benzene_variants_star_map)

    def test_get_item(self, benzene_variants_star_map):
        kc = KeyedChain.from_gufe(benzene_variants_star_map)

        assert kc[0] == kc._keyed_chain[0]
        assert kc[-1] == kc._keyed_chain[-1]
        assert kc[:] == kc._keyed_chain[:]


def test_datetime_to_json():
    d = datetime.datetime.fromisoformat('2023-05-05T09:06:43.699068')

    ser = json.dumps(d, cls=JSON_HANDLER.encoder)

    assert isinstance(ser, str)

    d2 = json.loads(ser, cls=JSON_HANDLER.decoder)

    assert isinstance(d2, datetime.datetime)
    assert d == d2

    reser = json.dumps(d2, cls=JSON_HANDLER.encoder)

    assert reser == ser
