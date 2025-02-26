from enum import Enum

import msgpack
import numpy as np

from gufe.tokenization import GufeTokenizable

class MPEXT(Enum):
    GUFETOKENIZABLE
    TIMESTAMP
    INT128
    SETTINGS
    UNIT
    NDARRAY
    PATH

def pack_default():
    raise NotImplementedError

def unpack_default():
    raise NotImplementedError

def packb(obj) -> bytes:
    return msgpack.packb(obj, default=default)

def unpackb(data: bytes):
    raise NotImplementedError
