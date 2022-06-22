from typing import NamedTuple, Dict

import hashlib

import importlib
import inspect


REMAPPED_CLASSES = {}
"""Mapping of old module/class names to new ones. Backward compatibility.

This should be a mapping of the Tuple[str, str] giving (old_module,
old_class) to Tuple[str, str] giving (new_module, new_class)
"""

class SerializationInfo(NamedTuple):
    """
    Parameters
    ----------
    bytes_data: byes
        bytes to be stored to disk for this object
    metadata: Dict
        metadata about this object, contains the following fields:

        * ``:path:``: The storage path for this object
        * ``:class:``: The qualname for this class of this object
        * ``:module:``: The module where the class is found
        * ``:md5:``: The md5 hexdigest for this object. In order to preserve
          this value across code versions, this hash is based on the
          arguments with non-default values, so it does not necessarily
          correspond to the md5 of the associated ``bytes_data``.
    """
    bytes_data: bytes
    metadata: Dict

    @property
    def path(self):
        return self.metadata[':path:']

    @property
    def md5(self):
        return self.metadata[':md5:']

    @staticmethod
    def create_metadata(obj, bytes_data, path_pattern):
        md5 = hashlib.md5(bytes_data).hexdigest()
        cls = obj.__class__
        return {
            ":path:": path_pattern.format(md5=md5[:10]),
            ":md5:": md5,
            ":class:": cls.__qualname__,
            ":module:": cls.__module__,
        }



def import_qualname(modname, qualname, remappings):
    if (modname, qualname) in remappings:
        modname, qualname = remappings[(modname, qualname)]

    result = importlib.import_module(modname)
    for name in qualname.split('.'):
        result = getattr(result, name)

    return result


def get_defaults_from_init(obj, dct):
    sig = inspect.signature(obj.__init__)
    defaultable = [param for param in sig.parameters.values()
                   if param.default is not inspect.Parameter.empty]
    defaults = [param.name for param in defaultable
                if param.default is dct[param.name]]
    return defaults

