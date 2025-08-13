# Parts of this file were originally part of OpenPathSampling/SimStore, which is
# distributed under the MIT license.
# Portions Copyright (c) 2014-2022 the contributors to OpenPathSampling
# Permissions are the same as those listed in the gufe LICENSE

# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import datetime
import functools
import json
import pathlib
from collections.abc import Callable, Iterable
from typing import Any, Dict, List, Optional, Tuple, Type, Union
from uuid import UUID

import numpy as np
from openff.units import DEFAULT_UNIT_REGISTRY

import gufe
from gufe.compression import zst_compress, zst_decompress
from gufe.settings.models import SettingsBaseModel


class JSONCodec:
    """Custom JSON encoding and decoding for non-default types.

    Parameters
    ----------
    cls : class
        Class for this codec. Assumes that all subclasses should be treated
        the same way. Can be ``None`` if ``is_my_obj`` and ``is_my_dict``
        are given.
    to_dict : Callable
        method that converts the object to a dictionary
    from_dict : Callable
        method that restores the object based on the dictionary made by
        to_dict
    is_my_obj : Optional[Callable]
        Method to determine whether the input object should be treated by
        this encoder. Default behavior is to use ``isinstance(cls)``, and to
        create a dict that also includes the class name and the module of
        the object.
    is_my_dict : Optional[Callable]
        Method to determine whether the input dictionary should be treated
        by this decoder. Default behavior assumes usage of the default
        ``is_my_obj``.
    """

    def __init__(
        self,
        cls: type | None,
        to_dict: Callable[[Any], dict],
        from_dict: Callable[[dict], Any],
        is_my_obj: Callable[[Any], bool] | None = None,
        is_my_dict=None,
    ):
        if is_my_obj is None:
            is_my_obj = self._is_my_obj

        if is_my_dict is None:
            is_my_dict = self._is_my_dict

        self.cls = cls
        self.to_dict = to_dict
        self.from_dict = from_dict
        self.is_my_obj = is_my_obj
        self.is_my_dict = is_my_dict

    def _is_my_dict(self, dct: dict) -> bool:
        expected = ["__class__", "__module__", ":is_custom:"]
        is_custom = all(exp in dct for exp in expected)
        return (
            is_custom
            and dct["__class__"] == self.cls.__name__  # type: ignore [union-attr]
            and dct["__module__"] == self.cls.__module__
        )

    def _is_my_obj(self, obj: Any) -> bool:
        return isinstance(obj, self.cls)  # type: ignore [arg-type]

    def default(self, obj: Any) -> Any:
        if self.is_my_obj(obj):
            dct = {}
            if self.cls:
                dct.update(
                    {
                        "__class__": obj.__class__.__qualname__,
                        "__module__": obj.__class__.__module__,
                        ":is_custom:": True,
                    }
                )
            # we let the object override __class__ and __module__ if needed
            dct.update(self.to_dict(obj))
            return dct
        return obj

    def object_hook(self, dct: dict) -> Any:
        if self.is_my_dict(dct):
            obj = self.from_dict(dct)
            return obj
        return dct


def custom_json_factory(
    coding_methods: Iterable[JSONCodec],
) -> tuple[type[json.JSONEncoder], type[json.JSONDecoder]]:
    """Create JSONEncoder/JSONDecoder for special types.

    Factory method. Dynamically creates classes that enable all the provided
    ``coding_methods``. Returns classes, not instances, as classes are used
    by the ``cls`` argument in ``json.loads`` / ``json.dumps``.

    Parameters
    ----------
    coding_methods : Iterable[JSONCodec]
        codecs to use

    Returns
    -------
    Tuple[Type[JSONEncoder], Type[JSONDecoder]]
        subclasses of JSONEncoder/JSONDecoder that use support the provided
        codecs
    """

    class CustomJSONEncoder(json.JSONEncoder):
        def default(self, obj):
            for coding_method in coding_methods:
                # If the coding method cannot handle this object, it returns
                # the object (unchanged). So if the object is changed, we
                # return that.
                result = coding_method.default(obj)
                if result is not obj:
                    return result

            # if none of our methods are useful, use standard approaches
            # (including providing standard error)
            return json.JSONEncoder.default(self, obj)

    class CustomJSONDecoder(json.JSONDecoder):
        def __init__(self, *args, **kwargs):
            # technically, JSONDecoder doesn't come with an object_hook
            # method, which is why we pass it to super here
            super().__init__(object_hook=self.object_hook, *args, **kwargs)

        def object_hook(self, dct):
            for coding_method in coding_methods:
                # If the coding method cannot handle this dict, it returns
                # the dict (unchanged). So if the dict is changed, we return
                # that.
                result = coding_method.object_hook(dct)
                if result is not dct:
                    return result

            # if none of our methods are useful, just return the dict
            return dct

    return (CustomJSONEncoder, CustomJSONDecoder)


class JSONSerializerDeserializer:
    r"""
    Tools to serialize and deserialize objects as JSON.

    This wrapper object is necessary so that we can register new codecs
    after the original initialization.

    Attributes
    ----------
    encoder:
        subclass of ``JSONEncoder``; use as ``json.dumps(obj, cls=encoder)``
    decoder:
        subclass of ``JSONDecoder``; use as ``json.loads(string,
        cls=decoder)``

    Parameters
    ----------
    codecs : list of :class:`.JSONCodec`\s
        codecs supported
    """

    def __init__(self, codecs: Iterable[JSONCodec]):
        self.codecs: list[JSONCodec] = []
        for codec in codecs:
            self.add_codec(codec)

        self.encoder, self.decoder = self._set_serialization()

    def _set_serialization(
        self,
    ) -> tuple[type[json.JSONEncoder], type[json.JSONDecoder]]:
        encoder, decoder = custom_json_factory(self.codecs)
        self._serializer = functools.partial(json.dumps, cls=encoder)
        self._deserializer = functools.partial(json.loads, cls=decoder)
        return encoder, decoder

    def add_codec(self, codec: JSONCodec):
        """Add a new codec to the supported codecs

        Parameters
        ----------
        codec : :class:`.JSONCodec`
            codec to add
        """
        if codec in self.codecs:
            return

        if codec is not None:
            self.codecs.append(codec)

        self.encoder, self.decoder = self._set_serialization()

    def serializer(self, obj: Any) -> str:
        """Callable that dumps to JSON"""
        return self._serializer(obj)

    def deserializer(self, string: str) -> Any:
        """Callable to loads JSON"""
        return self._deserializer(string)


def default_from_dict(dct):
    dct = dict(dct)  # make a copy
    module = dct.pop("__module__")
    qualname = dct.pop("__class__")
    del dct[":is_custom:"]

    cls = gufe.tokenization.get_class(module, qualname)
    return cls(**dct)


def inherited_is_my_dict(dct, cls):
    if not ("__module__" in dct and "__class__" in dct):
        return False
    dct = dict(dct)
    module = dct.pop("__module__")
    classname = dct.pop("__class__")
    stored = gufe.tokenization.get_class(module, classname)
    return cls in stored.mro()


def is_npy_dtype_dict(dct):
    expected = ["dtype", "bytes"]
    is_custom = all(exp in dct for exp in expected)
    return is_custom and ("shape" not in dct)


def is_openff_unit_dict(dct):
    expected = ["pint_unit_registry", "unit_name", ":is_custom:"]
    is_custom = all(exp in dct for exp in expected)
    return is_custom and dct["pint_unit_registry"] == "openff_units"


def is_openff_quantity_dict(dct):
    expected = ["pint_unit_registry", "magnitude", ":is_custom:", "unit"]
    is_custom = all(exp in dct for exp in expected)
    return is_custom and dct["pint_unit_registry"] == "openff_units"


def is_legacy_path_dict(dct: dict) -> bool:
    """This supports the case where python 3.12 needs to load python 3.13
    TODO: remove when 3.12 support is dropped
    """
    expected = ["__class__", "__module__", ":is_custom:", "path"]
    is_custom = all(exp in dct for exp in expected)
    return (
        is_custom
        and dct["__class__"] == "PosixPath"  # type: ignore [union-attr]
        and dct["__module__"] in ("pathlib", "pathlib._local")
    )


PATH_CODEC = JSONCodec(
    cls=pathlib.PosixPath,
    to_dict=lambda p: {"path": str(p)},
    from_dict=lambda dct: pathlib.PosixPath(dct["path"]),
    is_my_dict=is_legacy_path_dict,
)


BYTES_CODEC = JSONCodec(
    cls=bytes,
    to_dict=lambda obj: {"latin-1": zst_compress(obj).decode("latin-1")},
    from_dict=lambda dct: zst_decompress(dct["latin-1"].encode("latin-1")),
)


DATETIME_CODEC = JSONCodec(
    cls=datetime.datetime,
    to_dict=lambda obj: {"isotime": obj.isoformat()},
    from_dict=lambda dct: datetime.datetime.fromisoformat(dct["isotime"]),
)

# Note that this has inconsistent behaviour for some generic types
# which end up being handled by the default JSON encoder/decoder.
# The main example of this is np.float64 which will be turned into
# a float type on serialization.
NPY_DTYPE_CODEC = JSONCodec(
    cls=np.generic,
    to_dict=lambda obj: {
        "dtype": str(obj.dtype),
        "bytes": obj.tobytes(),
    },
    from_dict=lambda dct: np.frombuffer(dct["bytes"], dtype=np.dtype(dct["dtype"]))[0],
    is_my_obj=lambda obj: isinstance(obj, np.generic),
    is_my_dict=is_npy_dtype_dict,
)


NUMPY_CODEC = JSONCodec(
    cls=np.ndarray,
    to_dict=lambda obj: {
        "dtype": str(obj.dtype),
        "shape": list(obj.shape),
        "bytes": obj.tobytes(),
    },
    from_dict=lambda dct: np.frombuffer(dct["bytes"], dtype=np.dtype(dct["dtype"])).reshape(dct["shape"]),
)


SETTINGS_CODEC = JSONCodec(
    cls=SettingsBaseModel,
    to_dict=lambda obj: {field: getattr(obj, field) for field in obj.__fields__},
    from_dict=default_from_dict,
    is_my_dict=functools.partial(inherited_is_my_dict, cls=SettingsBaseModel),
)


OPENFF_QUANTITY_CODEC = JSONCodec(
    cls=None,
    to_dict=lambda obj: {
        "magnitude": obj.m,
        "unit": str(obj.u),
        ":is_custom:": True,
        "pint_unit_registry": "openff_units",
    },
    from_dict=lambda dct: dct["magnitude"] * DEFAULT_UNIT_REGISTRY.Quantity(dct["unit"]),
    is_my_obj=lambda obj: isinstance(obj, DEFAULT_UNIT_REGISTRY.Quantity),
    is_my_dict=is_openff_quantity_dict,
)


OPENFF_UNIT_CODEC = JSONCodec(
    cls=None,
    to_dict=lambda unit: {
        ":is_custom:": True,
        "pint_unit_registry": "openff_units",
        "unit_name": str(unit),
    },
    from_dict=lambda dct: DEFAULT_UNIT_REGISTRY(dct["unit_name"]).u,
    is_my_obj=lambda obj: isinstance(obj, DEFAULT_UNIT_REGISTRY.Unit),
    is_my_dict=is_openff_unit_dict,
)


UUID_CODEC = JSONCodec(
    cls=UUID,
    to_dict=lambda p: {"uuid": str(p)},
    from_dict=lambda dct: UUID(dct["uuid"]),
)
