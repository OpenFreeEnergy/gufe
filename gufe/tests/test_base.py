import pytest

from gufe.base import GufeTokenizable, GufeKey, tokenize


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


class TestGufeTokenizable:

    def test_init(self):
        leaf = Leaf("foo")
        bar = Leaf(leaf)

        cont = Container(bar, [leaf, 0], {"leaf": leaf, "a": "b"})


    def test_to_dict(self):
        ...

    def test_to_dict_roundtrip(self):
        ...

    def test_to_keyed_dict(self):
        ...

    def test_to_keyed_dict_roundtrip(self):
        ...

    def test_to_shallow_dict(self):
        ...

    def test_to_shallow_dict_roundtrip(self):
        ...
