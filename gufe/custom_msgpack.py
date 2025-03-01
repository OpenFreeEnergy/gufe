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
    LARGEINT = 0
    PATH = 1
    UNIT = 2
    NDARRAY = 3
    NPGENERIC = 4
    GUFETOKENIZABLE = 5
    SETTINGS = 6
    UUID = 7
    DATETIME = 8


def pack_largeint(large_int: int) -> bytes:
    return msgpack.packb(str(large_int))


def pack_default(obj) -> msgpack.ExtType:
    """For non-standard datatypes, create an ExtType containing packed
    data."""
    match obj:
        case GufeTokenizable():
            gt_payload: bytes = msgpack.packb(obj.to_keyed_chain(), default=pack_default)
            return msgpack.ExtType(MPEXT.GUFETOKENIZABLE, gt_payload)
        case datetime():
            dt_payload: bytes = msgpack.packb(obj.isoformat())
            return msgpack.ExtType(MPEXT.DATETIME, dt_payload)
        case SettingsBaseModel():
            settings_data = {field: getattr(obj, field) for field in obj.__fields__}
            settings_data.update({"__class__": obj.__class__.__qualname__, "__module__": obj.__class__.__module__})
            settings_payload: bytes = msgpack.packb(settings_data, default=pack_default)
            return msgpack.ExtType(MPEXT.SETTINGS, settings_payload)
        case pint.registry.UnitRegistry.Quantity():
            unit_data: list = [obj.m, str(obj.u), "openff_units"]
            unit_payload: bytes = msgpack.packb(unit_data, default=pack_default)
            return msgpack.ExtType(
                MPEXT.UNIT,
                unit_payload,
            )
        case Path():
            return msgpack.ExtType(MPEXT.PATH, msgpack.packb(str(obj)))
        case np.generic():
            npg_payload: bytes = msgpack.packb([str(obj.dtype), obj.tobytes()], default=pack_default)
            return msgpack.ExtType(MPEXT.NPGENERIC, npg_payload)
        case np.ndarray():
            npa_payload: bytes = msgpack.packb([str(obj.dtype), list(obj.shape), obj.tobytes()], default=pack_default)
            return msgpack.ExtType(MPEXT.NDARRAY, npa_payload)
        case int():
            try:
                return msgpack.packb(obj)
            except OverflowError:
                return msgpack.ExtType(MPEXT.LARGEINT, pack_largeint(obj))
        case uuid.UUID():
            # since a UUID is really a 128 bit int, rely on the LARGEINT extension type downstream
            uuid_payload: bytes = msgpack.packb(int(obj), default=pack_default)
            return msgpack.ExtType(MPEXT.UUID, uuid_payload)


def unpack_default(code: int, data: bytes):
    """For non-standard datatypes, unpack into the appropriate structures."""

    try:
        extension_type = MPEXT(code)
    except ValueError as e:
        raise ValueError(f"Found an unknown extension code: {code}")

    match extension_type:
        case MPEXT.GUFETOKENIZABLE:
            return gufe.tokenization.KeyedChain(
                msgpack.unpackb(data, ext_hook=unpack_default, strict_map_key=False)
            ).to_gufe()
        case MPEXT.SETTINGS:
            settings_data = msgpack.unpackb(data, ext_hook=unpack_default)
            module = settings_data.pop("__module__")
            qualname = settings_data.pop("__class__")
            cls = gufe.tokenization.get_class(module, qualname)
            return cls(**settings_data)
        case MPEXT.UNIT:
            magnitude, unit, _ = msgpack.unpackb(data, ext_hook=unpack_default)
            return magnitude * DEFAULT_UNIT_REGISTRY.Quantity(unit)
        case MPEXT.NDARRAY:
            _dtype, _shape, _bytes = msgpack.unpackb(data)
            return np.frombuffer(_bytes, dtype=np.dtype(_dtype)).reshape(_shape)
        case MPEXT.LARGEINT:
            string_rep_data = msgpack.unpackb(data)
            return int(string_rep_data)
        case MPEXT.PATH:
            string_rep_data = msgpack.unpackb(data)
            return Path(string_rep_data)
        case MPEXT.DATETIME:
            isoformat = msgpack.unpackb(data)
            return datetime.fromisoformat(isoformat)
        case MPEXT.UUID:
            uuid_data = msgpack.unpackb(data, ext_hook=unpack_default)
            return uuid.UUID(int=int(uuid_data))
        case MPEXT.NPGENERIC:
            _dtype, _bytes = msgpack.unpackb(data, ext_hook=unpack_default)
            return np.frombuffer(_bytes, dtype=np.dtype(_dtype))


def packb(obj) -> bytes:
    return msgpack.packb(obj, default=pack_default)


def unpackb(data: bytes):
    return msgpack.unpackb(data, ext_hook=unpack_default)
