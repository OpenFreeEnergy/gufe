import abc
import uuid
import sys
import threading
import hashlib
import inspect
from packaging.version import parse as parse_version
from typing import Dict, Any, Callable


class GufeTokenizableMixin:
    """Mixin for all tokenizble tokenizeable gufe objects.

    """

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
        return self._key

    @property
    def defaults(self):
        sig = inspect.signature(self.__init__)

        defaults = {
            param.name: param.default for param in sig.parameters.values()
            if param.default is not inspect.Parameter.empty
        }

        return defaults

    def _keyencode_dependencies(self, d):
        for key, value in d.items():
            if isinstance(value, dict):
                self._keyencode_dependencies(value)
            elif isinstance(value, list):
                for i, item in enumerate(value):
                    if hasattr(item, '_gufe_tokenize'):
                        d[key][i] = item
            else:
                if hasattr(value, '_gufe_tokenize'):
                    d[key] = value.token

        return d


class GufeKey(str):
    def __repr__(self):
        return f"<GufeKey('{str(self)}')>"


STUBBED_OBJECT_REGISTRY: Dict[str, GufeTokenizableMixin] = {}
"""Registry of token-stubbed objects.

Used to avoid duplication of tokenizable `gufe` objects in memory when deserialized.
Each key is a token, each value the token-stubbed version of the object
corresponding to that token.

"""

#TODO see if you can proceed without this first
#STUBBED_OBJECT_REGISTRY_REVERSE: Dict[GufeTokenizableMixin, str] = {}
#"""Reverse registry of token-stubbed objects.
#
#Used to Avoid duplication of tokenizable `gufe` objects in memory when deserialized.
#Each key is a token-stubbed version of a `gufe` object, and each value is a 
#token corresponding to that object.
#
#"""


# FROM dask.base

# Pass `usedforsecurity=False` for Python 3.9+ to support FIPS builds of Python
_PY_VERSION = parse_version(".".join(map(str, sys.version_info[:3])))
_md5: Callable
if _PY_VERSION >= parse_version("3.9"):

    def _md5(x, _hashlib_md5=hashlib.md5):
        return _hashlib_md5(x, usedforsecurity=False)

else:
    _md5 = hashlib.md5

NORMALIZERS = {}

def normalize_object(o):
    method = getattr(o, "_gufe_tokenize", None)
    if method is not None:
        return method()
    else:
        raise ValueError("Cannot normalize without `_gufe_tokenize` method.")


def tokenize(arg):
    """Deterministic token"""
    hasher = _md5(str(normalize_object(arg)).encode())
    return hasher.hexdigest()


