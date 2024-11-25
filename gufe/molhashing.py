# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import io

import numpy as np


def serialize_numpy(arr: np.ndarray) -> str:
    """Convert a numpy array to a string for serialisation

    SeeAlso
    -------
    :func:`.deserialize_numpy`
    """
    npbytes = io.BytesIO()
    np.save(npbytes, arr, allow_pickle=False)
    npbytes.seek(0)
    # latin-1 or base64? latin-1 is fewer bytes, but arguably worse on eyes
    return npbytes.read().decode("latin-1")


def deserialize_numpy(arr_str: str) -> np.ndarray:
    """Convert a string representation to a numpy array

    SeeAlso
    -------
    :func:`.serialize_numpy`
    """
    npbytes = io.BytesIO(arr_str.encode("latin-1"))
    npbytes.seek(0)
    return np.load(npbytes)
