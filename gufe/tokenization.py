# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""
The machinery for tokenizing gufe objects live in this module.
"""
import abc
import hashlib
import importlib
import inspect
import json
import logging
import networkx as nx
import re
import warnings
import weakref
from itertools import chain
from typing import Any, Union, List, Tuple, Dict, Generator
from typing_extensions import Self

from gufe.custom_codecs import (
    BYTES_CODEC,
    DATETIME_CODEC,
    NPY_DTYPE_CODEC,
    NUMPY_CODEC,
    OPENFF_QUANTITY_CODEC,
    OPENFF_UNIT_CODEC,
    PATH_CODEC,
    SETTINGS_CODEC,
    UUID_CODEC,
)
from gufe.custom_json import JSONSerializerDeserializer

_default_json_codecs = [
    PATH_CODEC,
    NUMPY_CODEC,
    NPY_DTYPE_CODEC,
    BYTES_CODEC,
    DATETIME_CODEC,
    SETTINGS_CODEC,
    OPENFF_UNIT_CODEC,
    OPENFF_QUANTITY_CODEC,
    UUID_CODEC,
]
JSON_HANDLER = JSONSerializerDeserializer(_default_json_codecs)

# maps qualified name strings to the class
TOKENIZABLE_CLASS_REGISTRY: dict[tuple[str, str], "GufeTokenizable"] = {}
""":noindex:"""

# maps fully qualified names to a new location
# e.g. if we did a rename:
# ('gufe', 'ProteinComponent') -> ('gufe', 'BiopolymerComponent')
REMAPPED_CLASSES: dict[tuple[str, str], tuple[str, str]] = {}
""":noindex:"""


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
        key = instance.key
        TOKENIZABLE_REGISTRY.setdefault(key, instance)
        instance = TOKENIZABLE_REGISTRY[key]
        return instance

    def __init__(cls, clsname, bases, attrs):
        """
        Restore the signature of __init__ or __new__
        """
        if inspect.signature(cls.__new__) != inspect.signature(object.__new__):
            sig = inspect.signature(cls.__new__)
        elif inspect.signature(cls.__init__) != inspect.signature(object.__init__):
            sig = inspect.signature(cls.__init__)
        else:
            # No __new__ or __init__ method defined
            return super().__init__(clsname, bases, attrs)

        # Remove the first parameter (cls/self)
        parameters = tuple(sig.parameters.values())
        cls.__signature__ = sig.replace(parameters=parameters[1:])

        return super().__init__(clsname, bases, attrs)


class _ABCGufeClassMeta(_GufeTokenizableMeta, abc.ABCMeta):
    # required to make use of abc.ABC in classes that use _ComponentMeta
    # metaclass: see https://stackoverflow.com/questions/57349105/
    ...


class _GufeLoggerAdapter(logging.LoggerAdapter):
    """LoggerAdapter to insert the gufe key into contextual information.

    This allows logging users to use ``%(gufekey)s`` in their logging
    formatter strings, similarly to ``%(name)s`` or ``%(levelname)s``.

    For details on logger adapters, see the Python logging documentation:
    https://docs.python.org/3/library/logging.html#loggeradapter-objects

    Parameters
    ----------
    logger: :class:`logging.Logger`
        the logger for this class
    extra: :class:`.GufeTokenizable`
        the instance this adapter is associated with
    """
    def process(self, msg, kwargs):
        extra = kwargs.get('extra', {})
        if (extra_dict := getattr(self, '_extra_dict', None)) is None:
            try:
                gufekey = self.extra.key.split('-')[-1]
            except Exception:
                # no matter what happened, we have a bad key
                gufekey = "UNKNOWN"
                save_extra_dict = False
            else:
                save_extra_dict = True

            extra_dict = {
                'gufekey': gufekey
            }

            if save_extra_dict:
                self._extra_dict = extra_dict

        extra.update(extra_dict)
        kwargs['extra'] = extra
        return msg, kwargs


def new_key_added(dct, new_key, default):
    """Serialization migration: Add a new key to the dictionary.

    This can be used in when writing a serialization migration (see
    :meth:`GufeTokenizable.serialization_migration` ) where a new key has
    been added to the object's representation (e.g., a new parameter has
    been added). In order to be migratable, the new key must have an
    associated default value.

    Parameters
    ----------
    dct : dict
        dictionary based on the old serialization version
    new_key: str
        name of the new key
    default: Any
        default value for the new key

    Returns
    -------
    dict:
        input dictionary modified to add the new key
    """
    dct[new_key] = default
    return dct


def old_key_removed(dct, old_key, should_warn):
    """Serialization migration: Remove an old key from the dictionary.

    This can be used in when writing a serialization migration (see
    :meth:`GufeTokenizable.serialization_migration` ) where a key has been
    removed from the object's serialized representation (e.g., an old
    parameter is no longer allowed). If a parameter has been removed, it is
    likely that you will want to warn the user that the parameter is no
    longer used: the ``should_warn`` option allows that.

    Parameters
    ----------
    dct : dict
        dictionary based on the old serialization version
    old_key : str
        name of the key that has been removed
    should_warn : bool
        whether to issue a warning for this (generally recommended)

    Returns
    -------
    dict:
        input dictionary modified to remove the old key
    """
    if should_warn:
        # TODO: this should be put elsewhere so that the warning can be more
        # meaningful (somewhere that knows what class we're recreating)
        warnings.warn(f"Outdated serialization: '{old_key}', with value "
                      f"'{dct[old_key]}' is no longer used in this object")

    del dct[old_key]
    return dct


def key_renamed(dct, old_name, new_name):
    """Serialization migration: Rename a key in the dictionary.

    This can be used in when writing a serialization migration (see
    :meth:`GufeTokenizable.serialization_migration` ) where a key has been
    renamed (e.g., a parameter name has changed).

    Parameters
    ----------
    dct : dict
        dictionary based on the old serialization version
    old_name : str
        name of the key in the old serialization representation
    new_name : str
        name of the key in the new serialization representation

    Returns
    -------
    dict:
        input dictionary modified to rename the key from the old name to the
        new one
    """
    dct[new_name] = dct.pop(old_name)
    return dct


def _label_to_parts(label):
    """Helper to split labels used for nested objects into parts.

    See :func:`.nested_key_moved` for a description of the label.
    """
    def _intify_if_possible(part):
        try:
            part = int(part)
        except ValueError:
            pass

        return part
    parts = [
        _intify_if_possible(p) for p in re.split('\.|\[|\]', label)
        if p != ""
    ]
    return parts


def _pop_nested(container, label):
    """Pop a nested object with the given label from a container.

    See :func:`.nested_key_moved` for a description of the label.
    """
    parts = _label_to_parts(label)
    current = container
    for part in parts[:-1]:
        current = current[part]

    return current.pop(parts[-1])


def _set_nested(container, label, value):
    """Set a nested object with the given label to the given value.

    See :func:`.nested_key_moved` for a description of the label.
    """
    parts = _label_to_parts(label)
    current = container
    for part in parts[:-1]:
        current = current[part]

    current[parts[-1]] = value


def nested_key_moved(dct, old_name, new_name):
    """Serialization migration: Move nested key to a new location.

    This can be used in when writing a serialization migration (see
    :meth:`GufeTokenizable.serialization_migration`) where a key that is
    nested in a structure of dicts/lists has been moved elsewhere. It uses
    labels that match Python namespace/list notations. That is, if ``dct``
    is the following dict::

        {'first': {'inner': ['list', 'of', 'words']}}

    then the label ``'first.inner[1]'`` would refer to the word ``'of'``.

    In that case, the following call::

        nested_key_moved(dct, 'first.inner[1]', 'second')

    would result in the dictionary::

        {'first': {'inner': ['list', 'words']}, 'second': 'of'}

    This is particular useful for things like protocol settings, which
    present as nested objects like this.

    Parameters
    ----------
    dct : dict
        dictionary based on the old serialization version
    old_name : str
        label for the old location (see above for description of label
        format)
    new_name : str
        label for the new location (see above for description of label
        format)

    Returns
    -------
    dict:
        input dictionary modified to move the value at the old location to
        the new location
    """
    val = _pop_nested(dct, old_name)
    _set_nested(dct, new_name, val)
    return dct


class GufeTokenizable(abc.ABC, metaclass=_ABCGufeClassMeta):
    """Base class for all tokenizeable gufe objects.

    Subclassing from this provides sorting, equality and hashing operators,
    provided that the class implements the `_to_dict` and `_from_dict` method.

    This extra work in serializing is important for hashes that are stable
    *across different Python sessions*.
    """
    @classmethod
    def _schema_version(cls) -> int:
        return 1

    def __repr__(self):
        return f"<{self.key}>"

    def __lt__(self, other):
        return self.key < other.key

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False

        return self.key == other.key

    def __hash__(self):
        return hash(self.key)

    def _gufe_tokenize(self):
        """Return a list of normalized inputs for `gufe.base.tokenize`.

        """
        return tokenize(self)
        # return normalize(self.to_keyed_dict(include_defaults=False))

    @classmethod
    def serialization_migration(cls, old_dict: dict, version: int) -> dict:
        """Migrate old serialization dicts to the current form.

        The input dict ``old_dict`` comes from some previous serialization
        version, given by ``version``. The output dict should be in the
        format of the current serialization dict.

        The recommended pattern to use looks like this:

        .. code::

            def serialization_migration(cls, old_dict, version):
                if version == 1:
                    ...  # do things for migrating version 1->2
                if version <= 2:
                    ...  # do things for migrating version 2->3
                if version <= 3:
                    ...  # do things for migrating version 3->4
                # etc

        This approach steps through each old serialization model on its way
        to the current version. It keeps code relatively minimal and
        readable.

        As a convenience, the following functions are available to simplify
        the various kinds of changes that are likely to occur in as
        serializtion versions change:

        * :func:`.new_key_added`
        * :func:`.old_key_removed`
        * :func:`.key_renamed`
        * :func:`.nested_key_moved`

        Parameters
        ----------
        old_dict : dict
            dict as received from a serialized form
        version: int
            the serialization version of ``old_dict``

        Returns
        -------
        dict :
            serialization dict suitable for the current implementation of
            ``from_dict``.

        """
        return old_dict

    @property
    def logger(self):
        """Return logger adapter for this instance"""
        if (adapter := getattr(self, '_logger', None)) is None:
            cls = self.__class__
            logname = "gufekey." + cls.__module__ + "." + cls.__qualname__
            logger = logging.getLogger(logname)
            adapter = _GufeLoggerAdapter(logger, self)
            self._logger = adapter
        return adapter

    @property
    def key(self):
        if not hasattr(self, '_key') or self._key is None:
            prefix = self.__class__.__qualname__
            token = self._gufe_tokenize()
            self._key = GufeKey(f"{prefix}-{token}")

        return self._key

    def _set_key(self, key: str):
        """Manually set the key.

        This should only be done when reloading an object whose key is not
        deterministic.

        Parameters
        ----------
        key : str
            contents of the GufeKey for this object
        """
        if old_key := getattr(self, '_key', None):
            TOKENIZABLE_REGISTRY.pop(old_key)

        self._key = GufeKey(key)
        TOKENIZABLE_REGISTRY.setdefault(self.key, self)

    @classmethod
    def defaults(cls):
        """Dict of default key-value pairs for this `GufeTokenizable` object.

        These defaults are stripped from the dict form of this object produced
        with ``to_dict(include_defaults=False)`` where default values are present.

        """
        defaults = cls._defaults()
        defaults[':version:'] = cls._schema_version()
        return defaults

    @classmethod
    @abc.abstractmethod
    def _defaults(cls):
        """This method should be overridden to provide the dict of defaults
        appropriate for the `GufeTokenizable` subclass.

        """
        sig = inspect.signature(cls.__init__)

        defaults = {
            param.name: param.default for param in sig.parameters.values()
            if param.default is not inspect.Parameter.empty
        }

        return defaults

    @abc.abstractmethod
    def _to_dict(self) -> dict:
        """This method should be overridden to provide the dict form of the
        `GufeTokenizable` subclass.

        `GufeTokenizable` instances should *not* be used as keys in dicts
        within this object; even though they are hashable, this makes
        serialization into e.g. JSON difficult.

        """
        raise NotImplementedError()

    @classmethod
    @abc.abstractmethod
    def _from_dict(cls, dct: dict):
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
            for key, value in self.defaults().items():
                if dct.get(key) == value:
                    dct.pop(key)

        return dct

    @classmethod
    def from_dict(cls, dct: dict):
        """Generate an instance from full dict representation.

        Parameters
        ----------
        dct : Dict
            A dictionary produced by `to_dict` to instantiate from.
            If an identical instance already exists in memory, it will be
            returned.  Otherwise, a new instance will be returned.

        """
        return dict_decode_dependencies(dct)

    def to_keyed_dict(self, include_defaults=True) -> dict:
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
            for key, value in self.defaults().items():
                if dct.get(key) == value:
                    dct.pop(key)

        return dct

    @classmethod
    def from_keyed_dict(cls, dct: dict):
        """Generate an instance from keyed dict representation.

        Parameters
        ----------
        dct : Dict
            A dictionary produced by `to_keyed_dict` to instantiate from.
            If an identical instance already exists in memory, it will be
            returned.  Otherwise, a new instance will be returned.

        """
        return key_decode_dependencies(dct)

    def to_shallow_dict(self) -> dict:
        """Generate shallow dict representation, with all referenced
        `GufeTokenizable` objects left intact.

        See also
        --------
        :meth:`GufeTokenizable.to_dict`
        :meth:`GufeTokenizable.to_keyed_dict`

        """
        return to_dict(self)

    @classmethod
    def from_shallow_dict(cls, dct: dict):
        """Generate an instance from shallow dict representation.

        Parameters
        ----------
        dct : Dict
            A dictionary produced by `to_shallow_dict` to instantiate from.
            If an identical instance already exists in memory, it will be
            returned.  Otherwise, a new instance will be returned.

        """
        return from_dict(dct)

    def copy_with_replacements(self, **replacements):
        """Make a modified copy of this object.

        Since GufeTokenizables are immutable, this is essentially a shortcut
        to mutate the object. Note that the keyword arguments it takes are
        based on keys of the dictionaries used in the the
        ``_to_dict``/``_from_dict`` cycle for this object; in most cases
        that is the same as parameters to ``__init__``, but not always.

        This will always return a *new* object in memory. So using
        ``obj.copy_with_replacements()`` (with no keyword arguments) is a
        way to create a shallow copy: the object is different in memory, but
        its attributes will be the same objects in memory as the original.

        Parameters
        ----------
        replacements: Dict
            keyword arguments with keys taken from the keys given by the
            output of this object's ``to_dict`` method.
        """
        dct = self._to_dict()
        if invalid := set(replacements) - set(dct):
            raise TypeError(f"Invalid replacement keys: {invalid}. "
                            f"Allowed keys are: {set(dct)}")

        dct.update(replacements)
        return self._from_dict(dct)


class GufeKey(str):
    def __repr__(self):   # pragma: no cover
        return f"<GufeKey('{str(self)}')>"

    def to_dict(self):
        return {':gufe-key:': str(self)}

    @property
    def prefix(self) -> str:
        """Commonly indicates a classname"""
        return self.split('-')[0]

    @property
    def token(self) -> str:
        """Unique hash of this key, typically a md5 value"""
        return self.split('-')[1]


def gufe_objects_from_shallow_dict(
    obj: Union[List, Dict, GufeTokenizable]
) -> List[GufeTokenizable]:
    """Find GufeTokenizables within a shallow dict.

    This function recursively looks through the list/dict structures encoding
    GufeTokenizables and returns list of all GufeTokenizables found
    within those structures, which may be potentially nested.

    Parameters
    ----------
    obj
        The input data structure to recursively traverse. For the initial call
        of this function, this should be the shallow dict of a GufeTokenizable.
        Input of a GufeTokenizable will immediately return a base case.

    Returns
    -------
    List[GufeTokenizable]
        All GufeTokenizables found in the shallow dict representation of a
        GufeTokenizable.

    """
    if isinstance(obj, GufeTokenizable):
        return [obj]

    elif isinstance(obj, list):
        return list(
            chain.from_iterable([gufe_objects_from_shallow_dict(item) for item in obj])
        )

    elif isinstance(obj, dict):
        return list(
            chain.from_iterable(
                [gufe_objects_from_shallow_dict(item) for item in obj.values()]
            )
        )

    return []


def gufe_to_digraph(gufe_obj):
    """Recursively construct a DiGraph from a GufeTokenizable.

    The DiGraph encodes the dependency structure of the GufeTokenizable on
    other GufeTokenizables.
    """
    graph = nx.DiGraph()
    shallow_dicts = {}

    def add_edges(o):
        # if we've made a shallow dict before, we've already added this one
        # and all its dependencies; return `None` to avoid going down the tree
        # again
        sd = shallow_dicts.get(o.key)
        if sd is not None:
            return None

        # if not, then we make the shallow dict only once, add it to our index,
        # add edges to dependencies, and return it so we continue down the tree
        sd = o.to_shallow_dict()

        shallow_dicts[o.key] = sd

        # add the object node in case there aren't any connections
        graph.add_node(o)
        connections = gufe_objects_from_shallow_dict(sd)

        for c in connections:
            graph.add_edge(o, c)

        return sd

    sd = add_edges(gufe_obj)
    _ = modify_dependencies(sd, add_edges, is_gufe_obj, mode="encode")

    return graph


class KeyedChain(object):
    """Keyed chain representation encoder of a GufeTokenizable.

    The keyed chain representation of a GufeTokenizable provides a
    topologically sorted list of gufe keys and GufeTokenizable keyed dicts
    that can be used to fully recreate a GufeTokenizable without the need for a
    populated TOKENIZATION_REGISTRY.

    The class wraps around a list of tuples containing the gufe key and the
    keyed dict form of the GufeTokenizable.

    Examples
    --------
    We can create a keyed chain representation from any GufeTokenizable, such
    as:

    >>> from gufe.tokenization import KeyedChain
    >>> s = SolventComponent()
    >>> keyed_chain = KeyedChain.gufe_to_keyed_chain_rep(s)
    >>> keyed_chain
    [('SolventComponent-26b4034ad9dbd9f908dfc298ea8d449f',
      {'smiles': 'O',
       'positive_ion': 'Na+',
       'negative_ion': 'Cl-',
       'ion_concentration': '0.15 molar',
       'neutralize': True,
       '__qualname__': 'SolventComponent',
       '__module__': 'gufe.components.solventcomponent',
       ':version:': 1})]

    And we can do the reverse operation as well to go from a keyed chain
    representation back to a GufeTokenizable:

    >>> KeyedChain(keyed_chain).to_gufe()
    SolventComponent(name=O, Na+, Cl-)

    """

    def __init__(self, keyed_chain):
        self._keyed_chain = keyed_chain

    @classmethod
    def from_gufe(cls, gufe_object: GufeTokenizable) -> Self:
        """Initialize a KeyedChain from a GufeTokenizable."""
        return cls(cls.gufe_to_keyed_chain_rep(gufe_object))

    def to_gufe(self) -> GufeTokenizable:
        """Initialize a GufeTokenizable."""
        gts: Dict[str, GufeTokenizable] = {}
        for gufe_key, keyed_dict in self:
            gt = key_decode_dependencies(keyed_dict, registry=gts)
            gts[gufe_key] = gt
        return gt

    @classmethod
    def from_keyed_chain_rep(cls, keyed_chain: List[Tuple[str, Dict]]) -> KeyedChain:
        """Initialize a KeyedChain from a keyed chain representation."""
        return cls(keyed_chain)

    def to_keyed_chain_rep(self) -> List[Tuple[str, Dict]]:
        """Return the keyed chain representation of this object."""
        return list(self)

    @staticmethod
    def gufe_to_keyed_chain_rep(
        gufe_object: GufeTokenizable,
    ) -> List[Tuple[str, Dict]]:
        """Create the keyed chain representation of a GufeTokenizable.

        This represents the GufeTokenizable as a list of two-element tuples
        containing, as their first and second elements, the gufe key and keyed
        dict form of the GufeTokenizable, respectively, and provides the
        underlying structure used in the KeyedChain class.

        Parameters
        ----------
        gufe_object
            The GufeTokenizable for which the KeyedChain is generated.

        Returns
        -------
        key_and_keyed_dicts
            The keyed chain representation of a GufeTokenizable.

        """
        key_and_keyed_dicts = [
            (str(gt.key), gt.to_keyed_dict())
            for gt in nx.topological_sort(gufe_to_digraph(gufe_object))
        ][::-1]
        return key_and_keyed_dicts

    def gufe_keys(self) -> Generator[str, None, None]:
        """Create a generator that iterates over the gufe keys in the KeyedChain."""
        for key, _ in self:
            yield key

    def keyed_dicts(self) -> Generator[Dict, None, None]:
        """Create a generator that iterates over the keyed dicts in the KeyedChain."""
        for _, _dict in self:
            yield _dict

    def __len__(self):
        return len(self._keyed_chain)

    def __iter__(self):
        return self._keyed_chain.__iter__()

    def __getitem__(self, index):
        return self._keyed_chain[index]


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


def modify_dependencies(obj: Union[dict, list], modifier, is_mine, mode, top=True):
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
        `GufeTokenizable` s or decode them from dict or key-encoded forms.
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
def to_dict(obj: GufeTokenizable) -> dict:
    dct = obj._to_dict()
    dct.update(module_qualname(obj))
    dct[':version:'] = obj._schema_version()
    return dct


def dict_encode_dependencies(obj: GufeTokenizable) -> dict:
    return modify_dependencies(
        obj.to_shallow_dict(),
        to_dict,
        is_gufe_obj,
        mode='encode',
        top=True
    )


def key_encode_dependencies(obj: GufeTokenizable) -> dict:
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

    if thing is None:  # -no-cov-
        return obj
    else:
        return thing


def _from_dict(dct: dict) -> GufeTokenizable:
    module = dct.pop('__module__')
    qualname = dct.pop('__qualname__')
    version = dct.pop(':version:', 1)

    cls = get_class(module, qualname)
    dct = cls.serialization_migration(dct, version)
    return cls._from_dict(dct)


def dict_decode_dependencies(dct: dict) -> GufeTokenizable:
    return from_dict(
        modify_dependencies(dct, from_dict, is_gufe_dict, mode='decode', top=True)
    )


def key_decode_dependencies(
    dct: dict,
    registry=TOKENIZABLE_REGISTRY
) -> GufeTokenizable:
    # this version requires that all dependent objects are already registered
    # responsibility of the storage system that uses this to do so
    dct = modify_dependencies(
        dct,
        lambda d: registry[GufeKey(d[":gufe-key:"])],
        is_gufe_key_dict,
        mode='decode',
        top=True
    )
    return from_dict(dct)


def get_all_gufe_objs(obj):
    """For GufeTokenizable obj, get all contained GufeTokenizables.

    This is useful when deduplicating GufeTokenizables for serialization.

    Parameters
    ----------
    obj: :class:`.GufeTokenizable`
        the container tokenizable

    Returns
    -------
    Set[GufeTokenizable]
        all contained GufeTokenizables
    """
    results = {obj}
    def modifier(o):
        results.add(o)
        return o.to_shallow_dict()

    _ = modify_dependencies(obj.to_shallow_dict(), modifier, is_gufe_obj,
                            mode='encode')
    return results


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
    # hasher = hashlib.md5(str(normalize(obj)).encode(), usedforsecurity=False)
    dumped = json.dumps(obj.to_keyed_dict(include_defaults=False),
                        sort_keys=True, cls=JSON_HANDLER.encoder)
    hasher = hashlib.md5(dumped.encode(), usedforsecurity=False)
    return hasher.hexdigest()
