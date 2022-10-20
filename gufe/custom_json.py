# This file was originally part of OpenPathSampling/SimStore, which is
# distributed under the MIT license.
# Portions Copyright (c) 2014-2022 the contributors to OpenPathSampling
# Permissions are the same as those listed in the gufe LICENSE

from typing import Tuple, Callable, Iterable, Dict, Any, List, Type
import json
import functools


class JSONCodec(object):
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
        cls: type,
        to_dict: Callable[[Any], Dict],
        from_dict: Callable[[Dict], Any],
        is_my_obj: Callable[[Any], bool] = None,
        is_my_dict=None
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
        expected = ['__class__', '__module__', ':is_custom:']
        is_custom = all(exp in dct for exp in expected)
        return (is_custom and dct['__class__'] == self.cls.__name__
                and dct['__module__'] == self.cls.__module__)

    def _is_my_obj(self, obj: Any) -> bool:
        return isinstance(obj, self.cls)

    def default(self, obj: Any) -> Any:
        if self.is_my_obj(obj):
            dct = {}
            if self.cls:
                dct.update({
                    '__class__': self.cls.__name__,
                    '__module__': self.cls.__module__,
                    ':is_custom:': True,
                })
            # we let the object override __class__ and __module__ if needed
            dct.update(self.to_dict(obj))
            return dct
        return obj

    def object_hook(self, dct: Dict) -> Any:
        if self.is_my_dict(dct):
            obj = self.from_dict(dct)
            return obj
        return dct


def custom_json_factory(
    coding_methods: Iterable[JSONCodec]
) -> Tuple[Type[json.JSONEncoder], Type[json.JSONDecoder]]:
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
            super(CustomJSONDecoder, self).__init__(
                object_hook=self.object_hook, *args, **kwargs
            )

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


class JSONSerializerDeserializer(object):
    """
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
    codecs : list of :class:`.JSONCodec`s
        codecs supported
    """
    def __init__(self, codecs: Iterable[JSONCodec]):
        self.codecs: List[JSONCodec] = []
        for codec in codecs:
            self.add_codec(codec)

        self.encoder, self.decoder = self._set_serialization()

    def _set_serialization(self) -> Tuple[Type[json.JSONEncoder],
                                          Type[json.JSONDecoder]]:
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
