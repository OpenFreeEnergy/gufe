import uuid
from datetime import datetime
from enum import IntEnum
from pathlib import Path

import msgpack
import numpy as np
import pint
from openff.units import DEFAULT_UNIT_REGISTRY

import gufe
from gufe.settings import SettingsBaseModel
from gufe.tokenization import GufeTokenizable


class MPEXT(IntEnum):
    TIMESTAMP = 0
    INT128 = 1
    PATH = 2
    UNIT = 3
    NDARRAY = 4
    NPGENERIC = 5
    GUFETOKENIZABLE = 6
    SETTINGS = 7
    UUID = 8


def pack_int128(large_int: int) -> bytes:
    return msgpack.packb(str(large_int))


def pack_default(obj) -> msgpack.ExtType:
    """For non-standard datatypes, create an ExtType containing packed
    data."""
    match obj:
        case GufeTokenizable():
            data = obj.to_keyed_chain()
            return msgpack.ExtType(MPEXT.GUFETOKENIZABLE.value, msgpack.packb(data, default=pack_default))
        case datetime():
            data = msgpack.ext.Timestamp.from_datetime(obj)
            return msgpack.ExtType(MPEXT.TIMESTAMP.value, data.to_bytes())
        case SettingsBaseModel():
            data = {field: getattr(obj, field) for field in obj.__fields__}
            data.update({"__class__": obj.__class__.__qualname__, "__module__": obj.__class__.__module__})
            return msgpack.ExtType(MPEXT.SETTINGS.value, msgpack.packb(data, default=pack_default))
        case pint.Quantity():

            return msgpack.ExtType(
                MPEXT.UNIT.value,
                msgpack.packb(
                    {
                        "magnitude": obj.m,
                        "unit": str(obj.u),
                        "pint_unit_registry": "openff_units",
                    },
                    default=pack_default,
                ),
            )
        case Path():
            return msgpack.ExtType(MPEXT.PATH.value, msgpack.packb(str(obj)))
        case np.generic():
            data = [str(obj.dtype), obj.tobytes()]
            return msgpack.ExtType(MPEXT.NPGENERIC.value, msgpack.packb(data, default=pack_default))
        case np.ndarray():
            # data = {"dtype": str(obj.dtype), "shape": list(obj.shape), "bytes": obj.tobytes()}
            data = [str(obj.dtype), list(obj.shape), obj.tobytes()]
            return msgpack.ExtType(MPEXT.NDARRAY.value, msgpack.packb(data))
        case int():
            try:
                return msgpack.packb(obj)
            except OverflowError:
                return msgpack.ExtType(MPEXT.INT128.value, pack_int128(obj))
        case uuid.UUID():
            # since a UUID is really a 128 bit int, rely on the INT128 extension type downstream
            data = msgpack.packb(int(obj), default=pack_default)
            return msgpack.ExtType(MPEXT.UUID.value, data)


def unpack_default(code: int, data: bytes):
    """For non-standard datatypes, unpack into the appropriate structures."""
    match MPEXT(code):
        case MPEXT.GUFETOKENIZABLE:
            return gufe.tokenization.KeyedChain(
                msgpack.unpackb(data, ext_hook=unpack_default, strict_map_key=False)
            ).to_gufe()
        case MPEXT.SETTINGS:
            data = msgpack.unpackb(data, ext_hook=unpack_default)
            try:
                module = data.pop("__module__")
                qualname = data.pop("__class__")
            except KeyError as e:
                raise e
            cls = gufe.tokenization.get_class(module, qualname)
            return cls(**data)
        case MPEXT.UNIT:
            data = msgpack.unpackb(data, ext_hook=unpack_default)
            unit_data = data["magnitude"] * DEFAULT_UNIT_REGISTRY.Quantity(data["unit"])
            return data["magnitude"] * DEFAULT_UNIT_REGISTRY.Quantity(data["unit"])
        case MPEXT.NDARRAY:
            _dtype, _shape, _bytes = msgpack.unpackb(data)
            return np.frombuffer(_bytes, dtype=np.dtype(_dtype)).reshape(_shape)
        case MPEXT.INT128:
            string_rep_data = msgpack.unpackb(data)
            return int(string_rep_data)
        case MPEXT.PATH:
            string_rep_data = msgpack.unpackb(data)
            return Path(string_rep_data)
        case MPEXT.TIMESTAMP:
            return msgpack.ext.Timestamp.from_bytes(data).to_datetime()
        case MPEXT.UUID:
            data = msgpack.unpackb(data, ext_hook=unpack_default)
            return uuid.UUID(int=int(data))
        case MPEXT.NPGENERIC:
            _dtype, _bytes = msgpack.unpackb(data, ext_hook=unpack_default)
            return np.frombuffer(_bytes, dtype=np.dtype(_dtype))
        case _:
            raise NotImplementedError(code)


def packb(obj) -> bytes:
    return msgpack.packb(obj, default=pack_default)


def unpackb(data: bytes):
    return msgpack.unpackb(data, ext_hook=unpack_default)
