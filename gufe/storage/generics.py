import json
from gufe.storage.utils import (
    import_qualname, SerializationInfo, REMAPPED_CLASSES
)


def _is_saved_metadata(dct):
    return (
        isinstance(dct, dict)
        and set(dct) == {":path:", ":md5:", ":module:", ":class:"}
    )

def generic_to_storage_ready(obj, path_pattern, defaults, *,
                             keys_to_replace=None, dict_rep_modifier=None,
                             key_to_attr=None):
    """
    Parameters
    ----------
    obj : Any
        object to serialize
    path_pattern : str
        string ready to be formatted; should include {md5}.
    default_keys : Iterable[str]
        keys of ``to_dict`` that have default values
    keys_to_replace : List[str]
        attribute names that should be replaced by the metadata for the
        object
    dict_rep_modifier: Optional[Callable[[Dict], Dict]]
        function to modify the default ``to_dict``, e.g., to add additional
        attributes; if not provided, no modification is made
    """
    if keys_to_replace is None:
        keys_to_replace = []

    # dict_rep will be the dict we turn into JSON
    dict_rep = obj.to_dict()

    if key_to_attr is None:
        key_to_attr = {k: k for k in dict_rep}

    # gather the attributes we need to replace and replace them
    dicts = [
        getattr(obj, keys_to_attr[key]).to_storage_ready()
        for key in keys_to_replace
    ]
    # ugly merge dicts; after py 3.10 use reduce(__or__, dicts)
    storage_ready = dict(sum([list(d.items()) for d in dicts], []))

    for key in keys_to_replace:
        attr_obj = getattr(obj, key_to_attr[key])
        dict_rep[key] = storage_ready[attr_obj].metadata

    # allow modifications to the dictionary (such as adding additional keys)
    # before writing out the bytes
    if dict_rep_modifier:
        dict_rep = dict_rep_modifier(dict_rep)

    bytes_data = json.dumps(dict_rep, sort_keys=True).encode("utf-8")

    dict_for_hash = {key: val for key, val in dict_rep.items()
                     if val != defaults[key]}

    hash_bytes = json.dumps(dict_for_hash, sort_keys=True).encode("utf-8")
    metadata = SerializationInfo.create_metadata(obj, hash_bytes,
                                                 path_pattern)

    # add this objects to the storage_ready dict and return it
    storage_ready[obj] = SerializationInfo(bytes_data, metadata)
    return storage_ready


def storage_bytes_to_dict(serialized_bytes, load_func):
    """
    Parameters
    ----------
    serialized_bytes : bytes
        serialization of this object; this is the ``bytes_data`` attribute
        of the :class:`.SerializationInfo`.
    load_func : Callable[[str], bytes]
        function to load a given path; this is typically
        ResultsClient.load_bytes
    """
    dct = json.loads(serialized_bytes.decode('utf-8'))
    for attr, val in dct.items():
        if _is_saved_metadata(val):
            val_cls = import_qualname(val[":module:"], val[":class:"],
                                      REMAPPED_CLASSES)
            val_bytes = load_func(val[":path:"])
            dct[attr] = val_cls.from_serialization_bytes(val_bytes)

    return dct
