# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

# custom_codecs.py: A place to keep various custom JSONCodec instances

import pathlib
import numpy as np

from gufe.custom_json import JSONCodec


PATH_CODEC = JSONCodec(
    cls=pathlib.Path,
    to_dict=lambda p: {'path': str(p)},
    from_dict=lambda dct: pathlib.Path(dct['path'])
)

BYTES_CODEC = JSONCodec(
    cls=bytes,
    to_dict=lambda obj: {'latin-1': obj.decode('latin-1')},
    from_dict=lambda dct: dct['latin-1'].encode('latin-1'),
)

NUMPY_CODEC = JSONCodec(
    cls=np.ndarray,
    to_dict=lambda obj: {
        'dtype': str(obj.dtype),
        'shape': list(obj.shape),
        'bytes': obj.tobytes()
    },
    from_dict=lambda dct: np.frombuffer(
        dct['bytes'], dtype=np.dtype(dct['dtype'])
    ).reshape(dct['shape'])
)


