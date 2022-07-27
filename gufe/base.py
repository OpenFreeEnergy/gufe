import abc
import uuid
import sys
import threading
import hashlib
import inspect
from copy import copy
from packaging.version import parse as parse_version
from typing import Dict, Any, Callable


TOKENIZABLE_CLASS_REGISTRY = {}

def register_tokenizable_class(cls):
    TOKENIZABLE_CLASS_REGISTRY[cls.__name__] = cls

def _hash(obj):
    return hash(obj.key)
    
class _GufeTokenizableMeta(type):
    def __new__(cls, name, bases, classdict):
        componentcls = super().__new__(cls, name, bases, classdict)
        register_tokenizable_class(componentcls)

        # since __hash__ isn't inheritable, monkey-patching is the best we can
        # do besides defining it in every subclass
        cls.__hash__ = _hash

        return componentcls


class _ABCGufeClassMeta(_GufeTokenizableMeta, abc.ABCMeta):
    # required to make use of abc.ABC in classes that use _ComponentMeta metaclass
    # see https://stackoverflow.com/questions/57349105/python-abc-inheritance-with-specified-metaclass 
    ...
        

class GufeTokenizable(metaclass=_ABCGufeClassMeta):
    """Base class for all tokenizeable gufe objects.

    """
    def __lt__(self, other):
        return hash(self) < hash(other)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented("Component comparisons must be of the same type")

        return self.to_dict() == other.to_dict()

    def __hash__(self):
        return hash(self.key)

    @property
    def token(self):
        if not hasattr(self, '_token') or self._token is None:
            self._token = tokenize(self)
        return self._token

    @property
    def key(self):
        if not hasattr(self, '_key') or self._key is None:
            prefix = type(self).__name__
            self._key = GufeKey(f"{prefix}-{self.token}")

        # add to registry if not already present
        TOKENIZABLE_REGISTRY[self._key] = self

        return self._key

    @property
    def defaults(self):
        sig = inspect.signature(self.__init__)

        defaults = {
            param.name: param.default for param in sig.parameters.values()
            if param.default is not inspect.Parameter.empty
        }

        return defaults

    # TODO: rewrite the below three functions to be fully recursive like 
    # _dictdecode_dependencies
    def _keyencode_dependencies(self, d):
        for key, value in d.items():
            if isinstance(value, dict):
                self._keyencode_dependencies(value)
            elif isinstance(value, list):
                for i, item in enumerate(value):
                    if hasattr(item, '_gufe_tokenize'):
                        d[key][i] = item.key
            else:
                if hasattr(value, '_gufe_tokenize'):
                    d[key] = value.key

        return d

    @classmethod
    def _keydecode_dependencies(cls, d):
        for key, value in d.items():
            if isinstance(value, dict):
                cls._keydecode_dependencies(value)
            elif isinstance(value, list):
                for i, item in enumerate(value):
                    if isinstance(item, str) and item in TOKENIZABLE_REGISTRY:
                        d[key][i] = TOKENIZABLE_REGISTRY[item]
            else:
                if isinstance(value, str) and value in TOKENIZABLE_REGISTRY:
                    d[key] = TOKENIZABLE_REGISTRY[value]

        return d

    def _dictencode_dependencies(self, d):
        for key, value in d.items():
            if isinstance(value, dict):
                self._dictencode_dependencies(value)
            elif isinstance(value, list):
                for i, item in enumerate(value):
                    if hasattr(item, 'to_dict'):
                        d[key][i] = item.to_dict()
            else:
                if hasattr(value, 'to_dict'):
                    d[key] = value.to_dict()

        return d

    @classmethod
    def _dictdecode_dependencies(cls, obj):
        if isinstance(obj, dict):
            if '__class__' in obj:
                obj = cls.from_dict(obj)
            else:
                for key, value in obj.items():
                    obj[key] = cls._dictdecode_dependencies(value)
        elif isinstance(obj, list):
            for i, item in enumerate(obj):
                obj[i] = cls._dictdecode_dependencies(item)

        return obj

    def _gufe_tokenize(self):
        """Return a list of normalized inputs for `gufe.base.tokenize`.

        """
        return sorted(self.to_dict(include_defaults=False).items(), key=str)

    def to_dict(self, include_defaults=True, keyencode_dependencies=False) -> dict:
        """Serialize to dict representation"""

        d = self._to_dict()
        d['__class__'] = self.__class__.__name__

        if keyencode_dependencies:
            d = self._keyencode_dependencies(d)
        else:
            d = self._dictencode_dependencies(d)

        if not include_defaults:
            for key, value in self.defaults.items():
                if d.get(key) == value:
                    d.pop(key)

        return d

    @classmethod
    def from_dict(cls, d: dict, keydecode_dependencies=False):
        """Deserialize from dict representation"""
        d = copy(d)
        component_type = d.pop('__class__', None)
        if component_type is None:
            raise KeyError("No `__class__` key found; unable to deserialize Component")

        if keydecode_dependencies:
            d = cls._keydecode_dependencies(d)
        else:
            d = cls._dictdecode_dependencies(d)

        obj = TOKENIZABLE_CLASS_REGISTRY[component_type]._from_dict(d)

        # if this object is already in memory, return existing object
        # new object will eventually be garbage collected
        if obj.key in TOKENIZABLE_REGISTRY:
            return TOKENIZABLE_REGISTRY[obj.key]
        else:
            return obj


class GufeKey(str):
    def __repr__(self):
        return f"<GufeKey('{str(self)}')>"



# TODO: may want to make this a weakref dict to avoid holding references here
TOKENIZABLE_REGISTRY: Dict[str, GufeTokenizable] = {}
"""Registry of tokenizable objects.

Used to avoid duplication of tokenizable `gufe` objects in memory when deserialized.
Each key is a token, each value the corresponding object.

"""


# FROM dask.base
# Pass `usedforsecurity=False` for Python 3.9+ to support FIPS builds of Python
_PY_VERSION = parse_version(".".join(map(str, sys.version_info[:3])))
_md5: Callable
if _PY_VERSION >= parse_version("3.9"):

    def _md5(x, _hashlib_md5=hashlib.md5):
        return _hashlib_md5(x, usedforsecurity=False)
else:
    _md5 = hashlib.md5


def normalize_object(o):
    try:
        o._gufe_tokenize()
    except AttributeError:
        raise ValueError("Cannot normalize object without `_gufe_tokenize` method.")


def tokenize(obj: GufeTokenizable):
    """Generate a deterministic, relatively-stable token from a `GufeTokenizable` object.

    Examples
    --------
    >>> from gufe import SolventComponent
    >>> s = SolventComponent()
    >>> tokenize(s)
    '6adf97f83acf6453d4a6a4b1070f3754'

    >>> tokenize(s) == tokenize(SolventComponent.from_dict(s.to_dict()))
    True

    """
    hasher = _md5(str(normalize_object(obj)).encode())
    return hasher.hexdigest()
