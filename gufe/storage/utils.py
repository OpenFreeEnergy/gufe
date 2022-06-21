from typing import NamedTuple, Dict

import hashlib

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
        * ``:md5:``: The md5 hexdigest for this object. In order to preserve
          this value across code versions, this hash is based on the
          arguments with non-default values, so it does not necessarily
          correspond to the md5 of the associated ``bytes_data``.
        * TODO: version or hash_attrs? something to clarify which inputs
          were to used to create the md5 metadata. I lean toward hash_attrs
          because it is more explicit (allows external things to easily
          recreate without needing to know gufe internals; should also make
          code easier to implement.) Requires that we add a _hash_attrs
          attribute to objects when we load them from storage.
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
        full_qualname = cls.__module__ + "." + cls.__qualname__
        return {
            ":path:": path_pattern.format(md5=md5),
            ":md5:": md5,
            ":class:": full_qualname,
        }


def generic_to_storage_ready(obj, path_pattern, *, attrs_to_replace=None,
                             dict_rep_modifier=None):
    """
    Parameters
    ----------
    obj : Any
        object to serialize
    path_pattern : str
        string ready to be formatted; should include {md5}.
    attrs_to_replace : List[str]
        attribute names that should be replaced by the metadata for the
        object
    dict_rep_modifier: Callable[[Dict], Dict]
        function to modify the default ``to_dict``, e.g., to add additional
        attributes
    """
    if attrs_to_replace is None:
        attrs_to_replace = []

    dicts = [
        getattr(obj, attr).to_storage_ready()
        for attr in attrs_to_replace
    ]
    storage_ready = dict_merge(*dicts)
    dict_rep = obj.to_dict()
    if dict_rep_modifier:
        dict_rep = dict_rep_modifier(dict_rep)

    bytes_data = _dict_to_bytes(dict_rep)
    storage_ready[obj] = SerializationInfo.from_bytes(obj, bytes_data,
                                                      path_pattern)
    metadata = _get_metadata(obj, byte_data, path_pattern)
    storage_ready[obj] = SerializationInfo(byte_data, metadata)
    return storage_ready

