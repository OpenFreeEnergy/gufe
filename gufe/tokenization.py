# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""
The machinery for tokenizing gufe objects live in this module.
"""
import abc
import datetime
import hashlib
import importlib
from pathlib import Path
import inspect
import copy
from typing import Dict, Any, Callable, Union, List, Tuple
import weakref
from gufe.custom_json import JSONSerializerDeserializer
from gufe.custom_codecs import PATH_CODEC, NUMPY_CODEC, BYTES_CODEC

_default_json_codecs = [PATH_CODEC, NUMPY_CODEC, BYTES_CODEC]
JSON_HANDLER = JSONSerializerDeserializer(_default_json_codecs)

# maps qualified name strings to the class
TOKENIZABLE_CLASS_REGISTRY: Dict[Tuple[str, str], "GufeTokenizable"] = {}
# maps fully qualified names to a new location
# e.g. if we did a rename:
# ('gufe', 'ProteinComponent') -> ('gufe', 'BiopolymerComponent')
REMAPPED_CLASSES: Dict[Tuple[str, str], Tuple[str, str]] = {}


def register_tokenizable_class(cls):
    TOKENIZABLE_CLASS_REGISTRY[(cls.__module__, cls.__qualname__)] = cls


class _GufeTokenizableMeta(type):
    def __new__(cls, name, bases, classdict):
        componentcls = super().__new__(cls, name, bases, classdict)
        register_tokenizable_class(componentcls)

        return componentcls

    def __call__(cls, *args, **kwargs):
        instance = super().__call__(*args, **kwargs)
        # add to registry if not already present
        TOKENIZABLE_REGISTRY.setdefault(instance.key, instance)
        return instance


class _ABCGufeClassMeta(_GufeTokenizableMeta, abc.ABCMeta):
    # required to make use of abc.ABC in classes that use _ComponentMeta
    # metaclass: see https://stackoverflow.com/questions/57349105/
    ...


class GufeTokenizable(abc.ABC, metaclass=_ABCGufeClassMeta):
    """Base class for all tokenizeable gufe objects.

    Subclassing from this provides sorting, equality and hashing operators,
    provided that the class implements the `_to_dict` and `_from_dict` method.

    This extra work in serializing is important for hashes that are stable
    *across different Python sessions*.
    """
    def __lt__(self, other):
        return self.key < other.key

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False

        return self.to_keyed_dict() == other.to_keyed_dict()

    def __hash__(self):
        return hash(self.key)

    def _gufe_tokenize(self):
        """Return a list of normalized inputs for `gufe.base.tokenize`.

        """
        return normalize(self.to_keyed_dict(include_defaults=False))

    @property
    def key(self):
        if not hasattr(self, '_key') or self._key is None:
            prefix = self.__class__.__qualname__
            token = tokenize(self)
            self._key = GufeKey(f"{prefix}-{token}")

        return self._key

    @property
    def defaults(self):
        """Dict of default key-value pairs for this `GufeTokenizable` object.

        These defaults are stripped from the dict form of this object produced
        with `to_dict(include_defaults=False) where default values are present.

        """
        return self._defaults()

    @abc.abstractmethod
    def _defaults(self):
        """This method should be overridden to provide the dict of defaults
        appropriate for the `GufeTokenizable` subclass.

        """
        sig = inspect.signature(self.__init__)

        defaults = {
            param.name: param.default for param in sig.parameters.values()
            if param.default is not inspect.Parameter.empty
        }

        return defaults

    @abc.abstractmethod
    def _to_dict(self) -> Dict:
        """This method should be overridden to provide the dict form of the
        `GufeTokenizable` subclass.

        `GufeTokenizable` instances should *not* be used as keys in dicts
        within this object; even though they are hashable, this makes
        serialization into e.g. JSON difficult.

        """
        raise NotImplementedError()

    @classmethod
    @abc.abstractmethod
    def _from_dict(cls, dct: Dict):
        """This method should be overridden to receive the dict form of the
        `GufeTokenizable` subclass and generate an instance from it.

        """
        raise NotImplementedError()

    def to_dict(self, include_defaults=True) -> dict:
        """Generate full dict representation, with all referenced
        `GufeTokenizable` objects also given in full dict representations.

        Parameters
        ----------
        include_defaults : bool
            If `False`, strip keys from dict representation with values equal
            to those in `defaults`.

        See also
        --------
        :meth:`GufeTokenizable.to_shallow_dict`
        :meth:`GufeTokenizable.to_keyed_dict`

        """
        dct = dict_encode_dependencies(self)

        if not include_defaults:
            for key, value in self.defaults.items():
                if dct.get(key) == value:
                    dct.pop(key)

        return dct

    @classmethod
    def from_dict(cls, dct: Dict):
        """Generate an instance from full dict representation.

        Parameters
        ----------
        dct : Dict
            A dictionary produced by `to_dict` to instantiate from.
            If an identical instance already exists in memory, it will be
            returned.  Otherwise, a new instance will be returned.

        """
        return dict_decode_dependencies(dct)

    def to_keyed_dict(self, include_defaults=True) -> Dict:
        """Generate keyed dict representation, with all referenced
        `GufeTokenizable` objects given in keyed representations.

        A keyed representation of an object is a dict of the form:

            {':gufe-key:': <GufeTokenizable.key>}

        These function as stubs to allow for serialization and storage of
        `GufeTokenizable` objects with minimal duplication.

        The original object can be re-assembled with `from_keyed_dict`.

        See also
        --------
        :meth:`GufeTokenizable.to_dict`
        :meth:`GufeTokenizable.to_shallow_dict`

        """
        dct = key_encode_dependencies(self)

        if not include_defaults:
            for key, value in self.defaults.items():
                if dct.get(key) == value:
                    dct.pop(key)

        return dct

    @classmethod
    def from_keyed_dict(cls, dct: Dict):
        """Generate an instance from keyed dict representation.

        Parameters
        ----------
        dct : Dict
            A dictionary produced by `to_keyed_dict` to instantiate from.
            If an identical instance already exists in memory, it will be
            returned.  Otherwise, a new instance will be returned.

        """
        return key_decode_dependencies(dct)

    def to_shallow_dict(self) -> Dict:
        """Generate shallow dict representation, with all referenced
        `GufeTokenizable` objects left intact.

        See also
        --------
        :meth:`GufeTokenizable.to_dict`
        :meth:`GufeTokenizable.to_keyed_dict`

        """
        return to_dict(self)

    @classmethod
    def from_shallow_dict(cls, dct: Dict):
        """Generate an instance from shallow dict representation.

        Parameters
        ----------
        dct : Dict
            A dictionary produced by `to_shallow_dict` to instantiate from.
            If an identical instance already exists in memory, it will be
            returned.  Otherwise, a new instance will be returned.

        """
        return from_dict(dct)


class GufeKey(str):
    def __repr__(self):   # pragma: no cover
        return f"<GufeKey('{str(self)}')>"

    def to_dict(self):
        return {':gufe-key:': str(self)}


# TOKENIZABLE_REGISTRY: Dict[str, weakref.ref[GufeTokenizable]] = {}
TOKENIZABLE_REGISTRY: weakref.WeakValueDictionary[str, GufeTokenizable] = weakref.WeakValueDictionary()
"""Registry of tokenizable objects.

Used to avoid duplication of tokenizable `gufe` objects in memory when
deserialized.  Each key is a token, each value the corresponding object.

We use a `weakref.WeakValueDictionary` here to avoid holding references to
objects that are no longer referenced anywhere else.

"""


def module_qualname(obj):
    return {'__qualname__': obj.__class__.__qualname__,
            '__module__': obj.__class__.__module__}


def is_gufe_obj(obj: Any):
    return isinstance(obj, GufeTokenizable)


def is_gufe_dict(dct: Any):
    return (isinstance(dct, dict) and '__qualname__' in dct
            and '__module__' in dct)


def is_gufe_key_dict(dct: Any):
    return isinstance(dct, dict) and ":gufe-key:" in dct


# conveniences to get a class from module/class name
def import_qualname(modname: str, qualname: str, remappings=REMAPPED_CLASSES):
    if (qualname is None) or (modname is None):
        raise ValueError("`__qualname__` or `__module__` cannot be None; "
                         f"unable to identify object {modname}.{qualname}")

    if (modname, qualname) in remappings:
        modname, qualname = remappings[(modname, qualname)]

    result = importlib.import_module(modname)
    for name in qualname.split('.'):
        result = getattr(result, name)

    return result


def get_class(module: str, qualname: str):
    key = module, qualname
    try:
        return TOKENIZABLE_CLASS_REGISTRY[key]
    except KeyError:
        cls = import_qualname(module, qualname, REMAPPED_CLASSES)
        TOKENIZABLE_CLASS_REGISTRY[key] = cls
        return cls


def modify_dependencies(obj: Union[Dict, List], modifier, is_mine, mode, top=True):
    """
    Parameters
    ----------
    obj : Dict or List
        Dictionary or list to traverse. Assumes that only mappings are dict and
        only iterables are list, and that no gufe objects are in the keys of
        dicts
    modifier : Callable[[GufeTokenizable], Any]
        Function that modifies any GufeTokenizable found
    is_mine : Callable[Any, bool]
        Function that determines whether the given object should be
        subjected to the modifier
    mode : {'encode', 'decode'}
        Whether this function is being used to encode a set of
        `GufeTokenizable`s or decode them from dict or key-encoded forms.
        Required to determine when to modify objects found in nested dict/list.
    top : bool
        If `True`, skip modifying `obj` itself; needed for recursive use to
        avoid early stopping on `obj`.
    """
    if is_mine(obj) and not top and mode == 'encode':
        obj = modifier(obj)

    if isinstance(obj, dict):
        obj = {key: modify_dependencies(value, modifier, is_mine, mode=mode, top=False)
               for key, value in obj.items()}

    elif isinstance(obj, list):
        obj = [modify_dependencies(item, modifier, is_mine, mode=mode, top=False)
               for item in obj]

    if is_mine(obj) and not top and mode == 'decode':
        obj = modifier(obj)

    return obj


# encode options
def to_dict(obj: GufeTokenizable) -> Dict:
    dct = obj._to_dict()
    dct.update(module_qualname(obj))
    return dct


def dict_encode_dependencies(obj: GufeTokenizable) -> Dict:
    return modify_dependencies(
        obj.to_shallow_dict(),
        to_dict,
        is_gufe_obj,
        mode='encode',
        top=True
    )


def key_encode_dependencies(obj: GufeTokenizable) -> Dict:
    return modify_dependencies(
        obj.to_shallow_dict(),
        lambda obj: obj.key.to_dict(),
        is_gufe_obj,
        mode='encode',
        top=True
    )


# decode options
def from_dict(dct) -> GufeTokenizable:
    obj = _from_dict(dct)
    # When __new__ is called to create ``obj``, it should be added to the
    # TOKENIZABLE_REGISTRY. However, there seems to be some case (race
    # condition?) where this doesn't happen, leading to a KeyError inside
    # the dictionary if we use []. (When you drop into PDB and run the same
    # line that gave the error, you get the object back.) With ``get``,
    # ``thing`` becomes None, which is also what it would be if the weakref
    # was to a deleted object.
    thing = TOKENIZABLE_REGISTRY.get(obj.key)

    if thing is None:
        return obj
    else:
        return thing


def _from_dict(dct: Dict) -> GufeTokenizable:
    module = dct.pop('__module__')
    qualname = dct.pop('__qualname__')

    cls = get_class(module, qualname)
    return cls._from_dict(dct)


def dict_decode_dependencies(dct: Dict) -> GufeTokenizable:
    return from_dict(
        modify_dependencies(dct, from_dict, is_gufe_dict, mode='decode', top=True)
    )


def key_decode_dependencies(dct: Dict) -> GufeTokenizable:
    # this version requires that all dependent objects are already registered
    # responsibility of the storage system that uses this to do so
    dct = modify_dependencies(
        dct,
        lambda d: TOKENIZABLE_REGISTRY[GufeKey(d[":gufe-key:"])],
        is_gufe_key_dict,
        mode='decode',
        top=True
    )
    return from_dict(dct)


## inspired by `dask.base`

# TODO: make this more efficient with a dispatch mechanism
# instead of a series of isinstance checks
def normalize(o):

    # dicts
    if isinstance(o, dict):
        dct = {key: normalize(value) 
               for key, value in o.items()}
        return normalize_dict(dct)

    # lists and tuples
    if isinstance(o, (list, tuple)):
        return list(map(normalize, o))

    # paths
    if isinstance(o, Path):
        return str(o)

    # GufeTokenizable
    method = getattr(o, "_gufe_tokenize", None)
    if method is not None:
        return method()

    # primitives we support
    if isinstance(o,
        (int,
        float,
        str,
        bytes,
        type(None),
        type,
        slice,
        complex,
        type(Ellipsis),
        datetime.date)):
            return o

    if isinstance(o, Exception):
        return 

    raise RuntimeError(
        f"Object {str(o)} cannot be deterministically hashed."
    )


def normalize_dict(d):
    return sorted(d.items(), key=str)


def tokenize(obj: GufeTokenizable) -> str:
    """Generate a deterministic, relatively-stable token from a
    `GufeTokenizable` object.

    Examples
    --------
    >>> from gufe import SolventComponent
    >>> s = SolventComponent()
    >>> tokenize(s)
    'e6eef7519854d35a5ce6c84136b3684c'

    >>> tokenize(s) == tokenize(SolventComponent.from_dict(s.to_dict()))
    True

    """
    hasher = hashlib.md5(str(normalize(obj)).encode(), usedforsecurity=False)
    return hasher.hexdigest()
