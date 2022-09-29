# This file was originally part of OpenPathSampling/SimStore, which is
# distributed under the MIT license.
# Portions Copyright (c) 2014-2022 the contributors to OpenPathSampling
# Permissions are the same as those listed in the gufe LICENSE

import json
import functools
import pathlib


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
    def __init__(self, codecs):
        self.codecs = []
        for codec in codecs:
            self.add_codec(codec)

        self.encoder, self.decoder = self._set_serialization()

    def _set_serialization(self):
        encoder, decoder = custom_json_factory(self.codecs)
        self._serializer = functools.partial(json.dumps, cls=encoder)
        self._deserializer = functools.partial(json.loads, cls=decoder)
        return encoder, decoder

    def add_codec(self, codec):
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


    def serializer(self, obj):
        """Callable that dumps to JSON"""
        return self._serializer(obj)

    def deserializer(self, string):
        """Callable to loads JSON"""
        return self._deserializer(string)


def custom_json_factory(coding_methods):
    """Create JSONEncoder/JSONDecoder for special types
    """
    class CustomJSONEncoder(json.JSONEncoder):
        def default(self, obj):
            for coding_method in coding_methods:
                result = coding_method.default(obj)
                if result is not obj:
                    return result
            return json.JSONEncoder.default(self, obj)

    class CustomJSONDecoder(json.JSONDecoder):
        def __init__(self, *args, **kwargs):
            super(CustomJSONDecoder, self).__init__(
                object_hook=self.object_hook, *args, **kwargs
            )

        def object_hook(self, dct):
            for coding_method in coding_methods:
                result = coding_method.object_hook(dct)
                if result is not dct:
                    return result
            return dct

    return (CustomJSONEncoder, CustomJSONDecoder)


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
    def __init__(self, cls, to_dict, from_dict, is_my_obj=None,
                 is_my_dict=None):
        if is_my_obj is None:
            is_my_obj = self._is_my_obj

        if is_my_dict is None:
            is_my_dict = self._is_my_dict

        self.cls = cls
        self.to_dict = to_dict
        self.from_dict = from_dict
        self.is_my_obj = is_my_obj
        self.is_my_dict = is_my_dict

    def _is_my_dict(self, dct):
        expected = ['__class__', '__module__', ':is_custom:']
        is_custom = all(exp in dct for exp in expected)
        if is_custom:
            return (dct['__class__'] == self.cls.__name__
                    and dct['__module__'] == self.cls.__module__)

    def _is_my_obj(self, obj):
        return isinstance(obj, self.cls)

    def default(self, obj):
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

    def object_hook(self, dct):
        if self.is_my_dict(dct):
            obj = self.from_dict(dct)
            return obj
        return dct


PATH_CODEC = JSONCodec(
    cls=pathlib.Path,
    to_dict=lambda p: {'path': str(p)},
    from_dict=lambda dct: pathlib.Path(dct['path'])
)
