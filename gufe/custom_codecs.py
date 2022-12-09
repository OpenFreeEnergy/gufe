# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

# custom_codecs.py: A place to keep various custom JSONCodec instances

import pathlib
import numpy as np

import gufe
from gufe.custom_json import JSONCodec
from gufe.settings.models import SettingsBaseModel
import openff.units
import functools

def default_from_dict(dct):
    dct = dict(dct)  # make a copy
    module = dct.pop('__module__')
    qualname = dct.pop('__class__')
    del dct[':is_custom:']

    cls = gufe.tokenization.get_class(module, qualname)
    return cls(**dct)


def inherited_is_my_dict(dct, cls):
    dct = dict(dct)
    module = dct.pop('__module__')
    classname = dct.pop('__class__')
    stored = gufe.tokenization.get_class(module, classname)
    return cls in stored.mro()


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

SETTINGS_CODEC = JSONCodec(
    cls=SettingsBaseModel,
    to_dict=lambda obj: {field: getattr(obj, field)
                         for field in obj.__fields__},
    from_dict=default_from_dict,
    is_my_dict=functools.partial(inherited_is_my_dict,
                                 cls=SettingsBaseModel)
)

OPENFF_QUANTITY_CODEC = JSONCodec(
    cls=openff.units.units.Quantity,
    to_dict=lambda obj: {'magnitude': obj.m, 'unit': obj.u},
    from_dict=lambda dct: dct['magnitude'] * dct['unit']
)

OPENFF_UNIT_CODEC = JSONCodec(
    cls=openff.units.units.Unit,
    to_dict=lambda obj: {'unit': str(obj)},
    from_dict=lambda dct: getattr(openff.units.unit, dct['unit']),
)
