# This file was originally part of OpenPathSampling/SimStore, which is
# distributed under the MIT license.
# Portions Copyright (c) 2014-2022 the contributors to OpenPathSampling
# Permissions are the same as those listed in the gufe LICENSE

import json
import pathlib

import numpy as np
import openff.units
from openff.units import unit
import pytest
from numpy import testing as npt
from uuid import uuid4
from gufe.custom_codecs import (
    BYTES_CODEC,
    NUMPY_CODEC,
    NPY_DTYPE_CODEC,
    OPENFF_QUANTITY_CODEC,
    OPENFF_UNIT_CODEC,
    PATH_CODEC,
    SETTINGS_CODEC,
    UUID_CODEC,
)
from gufe.custom_json import JSONSerializerDeserializer, custom_json_factory
from gufe import tokenization
from gufe.settings import models


class TestJSONSerializerDeserializer:
    def test_add_codec(self):
        # without bytes codec, can't serialize numpy
        serialization = JSONSerializerDeserializer([NUMPY_CODEC])
        obj = np.array([[1.0, 0.0], [2.0, 3.2]])
        with pytest.raises(TypeError):
            serialization.serializer(obj)
        # add the codec and it will work
        serialization.add_codec(BYTES_CODEC)
        serialized = serialization.serializer(obj)
        assert len(serialization.codecs) == 2
        reconstructed = serialization.deserializer(serialized)
        npt.assert_equal(obj, reconstructed)

    def test_add_existing_codec(self):
        serialization = JSONSerializerDeserializer([BYTES_CODEC])
        assert len(serialization.codecs) == 1
        serialization.add_codec(BYTES_CODEC)
        assert len(serialization.codecs) == 1


@pytest.mark.parametrize('obj', [
    np.array([[1.0, 0.0], [2.0, 3.2]]),
    np.float32(1.1)
])
@pytest.mark.parametrize('codecs', [
    [BYTES_CODEC, NUMPY_CODEC, NPY_DTYPE_CODEC],
    [NPY_DTYPE_CODEC, BYTES_CODEC, NUMPY_CODEC],
])
def test_numpy_codec_order_roundtrip(obj, codecs):
    serialization = JSONSerializerDeserializer(codecs)
    serialized = serialization.serializer(obj)
    reconstructed = serialization.deserializer(serialized)
    npt.assert_equal(obj, reconstructed)
    assert obj.dtype == reconstructed.dtype


class CustomJSONCodingTest:
    """Base class for testing codecs.

    In ``setup_method()``, user must define the following:

    * ``self.codec``: The codec to run
    * ``self.objs``: A list of objects to serialize
    * ``self.dcts``: A list of expected serilized forms of each object in
      ``self.objs``
    """

    def test_default(self):
        for (obj, dct) in zip(self.objs, self.dcts):
            assert self.codec.default(obj) == dct

    def test_object_hook(self):
        for (obj, dct) in zip(self.objs, self.dcts):
            assert self.codec.object_hook(dct) == obj

    def _test_round_trip(self, encoder, decoder):
        for (obj, dct) in zip(self.objs, self.dcts):
            json_str = json.dumps(obj, cls=encoder)
            reconstructed = json.loads(json_str, cls=decoder)
            assert reconstructed == obj
            json_str_2 = json.dumps(obj, cls=encoder)
            assert json_str == json_str_2

    def test_round_trip(self):
        encoder, decoder = custom_json_factory([self.codec])
        self._test_round_trip(encoder, decoder)

    def test_not_mine(self):
        # test that the default behavior is obeyed
        obj = {"test": 5}
        json_str = '{"test": 5}'
        encoder, decoder = custom_json_factory([self.codec])
        assert json.dumps(obj, cls=encoder) == json_str
        assert json.loads(json_str, cls=decoder) == obj


class TestNumpyCoding(CustomJSONCodingTest):
    def setup_method(self):
        self.codec = NUMPY_CODEC
        self.objs = [np.array([[1.0, 0.0], [2.0, 3.2]]), np.array([1, 0]),
                     np.array([1.0, 2.0, 3.0], dtype=np.float32)]
        shapes = [[2, 2], [2,], [3,]]
        dtypes = [str(arr.dtype) for arr in self.objs]  # may change by system?
        byte_reps = [arr.tobytes() for arr in self.objs]
        self.dcts = [
            {
                ":is_custom:": True,
                "__class__": "ndarray",
                "__module__": "numpy",
                "shape": shape,
                "dtype": dtype,
                "bytes": byte_rep,
            }
            for shape, dtype, byte_rep in zip(shapes, dtypes, byte_reps)
        ]

    def test_object_hook(self):
        # to get custom equality testing for numpy
        for (obj, dct) in zip(self.objs, self.dcts):
            reconstructed = self.codec.object_hook(dct)
            npt.assert_array_equal(reconstructed, obj)

    def test_round_trip(self):
        encoder, decoder = custom_json_factory([self.codec, BYTES_CODEC])
        for (obj, dct) in zip(self.objs, self.dcts):
            json_str = json.dumps(obj, cls=encoder)
            reconstructed = json.loads(json_str, cls=decoder)
            npt.assert_array_equal(reconstructed, obj)
            assert reconstructed.dtype == obj.dtype
            json_str_2 = json.dumps(obj, cls=encoder)
            assert json_str == json_str_2


class TestNumpyGenericCodec(TestNumpyCoding):
    def setup_method(self):
        self.codec = NPY_DTYPE_CODEC
        # Note that np.float64 is treated as a float by the
        # default json encode (and so returns a float not a numpy
        # object).
        self.objs = [np.bool_(True), np.float16(1.0), np.float32(1.0),
                     np.complex128(1.0),
                     np.clongdouble(1.0), np.uint64(1)]
        dtypes = [str(a.dtype) for a in self.objs]
        byte_reps = [a.tobytes() for a in self.objs]
        # Overly complicated extraction of the class name
        # to deal with the bool_ -> bool dtype class name problem
        classes = [str(a.__class__).split("'")[1].split('.')[1]
                   for a in self.objs]
        self.dcts = [
            {
                ":is_custom:": True,
                "__class__": classname,
                "__module__": "numpy",
                "dtype": dtype,
                "bytes": byte_rep,
            }
            for dtype, byte_rep, classname in zip(dtypes, byte_reps, classes)
        ]


class TestPathCodec(CustomJSONCodingTest):
    def setup_method(self):
        self.codec = PATH_CODEC
        self.objs = [
            pathlib.PosixPath("foo/bar"),
        ]
        self.dcts = [
            {
                ":is_custom:": True,
                "__class__": "PosixPath",
                "__module__": "pathlib",
                "path": "foo/bar",
            }
        ]


class TestSettingsCodec(CustomJSONCodingTest):
    def setup_method(self):
        self.codec = SETTINGS_CODEC
        self.objs = [
            models.Settings.get_defaults(),
        ]
        self.dcts = [
            {
                "__class__": "Settings",
                "__module__": "gufe.settings.models",
                ":is_custom:": True,
                "forcefield_settings": obj.forcefield_settings,
                "thermo_settings": obj.thermo_settings,
            }
            for obj in self.objs
        ]

        self.full_dump = [
            {
                "__class__": "Settings",
                "__module__": "gufe.settings.models",
                ":is_custom:": True,
                "forcefield_settings": {
                    "__class__": "OpenMMSystemGeneratorFFSettings",
                    "__module__": "gufe.settings.models",
                    ":is_custom:": True,
                    "constraints": "hbonds",
                    "rigid_water": True,
                    "hydrogen_mass": 3.0,
                    "forcefields": [
                        "amber/ff14SB.xml",
                        "amber/tip3p_standard.xml",
                        "amber/tip3p_HFE_multivalent.xml",
                        "amber/phosaa10.xml",
                    ],
                    "small_molecule_forcefield": "openff-2.0.0",
                    "nonbonded_method": "PME",
                    "nonbonded_cutoff": {':is_custom:': True,
                                         'magnitude': 1.0,
                                         'pint_unit_registry': 'openff_units',
                                         'unit': 'nanometer'},
                },
                "thermo_settings": {
                    "__class__": "ThermoSettings",
                    "__module__": "gufe.settings.models",
                    ":is_custom:": True,
                    "temperature": {":is_custom:": True,
                                    "magnitude": 300.0,
                                    "pint_unit_registry": "openff_units",
                                    "unit": "kelvin"},
                    "pressure": None,
                    "ph": None,
                    "redox_potential": None,
                },
            }
        ]
        self.required_codecs = [
            self.codec,
            OPENFF_QUANTITY_CODEC,
            OPENFF_UNIT_CODEC,
        ]

    def test_round_trip(self):
        encoder, decoder = custom_json_factory(self.required_codecs)
        self._test_round_trip(encoder, decoder)

    def test_full_dump(self):
        encoder, _ = custom_json_factory(self.required_codecs)
        for obj, dct in zip(self.objs, self.full_dump):
            as_str = json.dumps(obj, cls=encoder)
            as_dct = json.loads(as_str)  # turn off decoder here!
            assert dct == as_dct


class TestOpenFFQuantityCodec(CustomJSONCodingTest):
    def setup_method(self):
        self.codec = OPENFF_QUANTITY_CODEC
        self.objs = [
            openff.units.DEFAULT_UNIT_REGISTRY("1.0 * kg meter per second squared"),
        ]
        self.dcts = [
            {
                ":is_custom:": True,
                "magnitude": 1.0,
                "pint_unit_registry": "openff_units",
                "unit": "kilogram * meter / second ** 2",
            },
        ]


def test_openff_quantity_array_roundtrip():
    thing = unit.Quantity.from_list([
        (i + 1.0)*unit.kelvin for i in range(10)
    ])

    dumped = json.dumps(thing, cls=tokenization.JSON_HANDLER.encoder)

    returned = json.loads(dumped, cls=tokenization.JSON_HANDLER.decoder)

    assert returned.u == thing.u
    assert (returned.m == thing.m).all()


class TestOpenFFUnitCodec(CustomJSONCodingTest):
    def setup_method(self):
        self.codec = OPENFF_UNIT_CODEC
        self.objs = [
            openff.units.unit.amu,
        ]
        self.dcts = [
            {
                ":is_custom:": True,
                "pint_unit_registry": "openff_units",
                "unit_name": "unified_atomic_mass_unit",
            }
        ]

class TestUUIDCodec(CustomJSONCodingTest):
    def setup_method(self):
        self.codec = UUID_CODEC
        self.objs = [
            uuid4()
        ]
        self.dcts = [
            {
                ":is_custom:": True,
                "__class__": "UUID",
                "__module__": "uuid",
                "uuid": f"{str(self.objs[0])}",
            }
        ]