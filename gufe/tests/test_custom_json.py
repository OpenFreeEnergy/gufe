# This file was originally part of OpenPathSampling/SimStore, which is
# distributed under the MIT license.
# Portions Copyright (c) 2014-2022 the contributors to OpenPathSampling
# Permissions are the same as those listed in the gufe LICENSE

import json
import pathlib

import numpy as np
import openff.units
import pytest
from numpy import testing as npt

from gufe.custom_codecs import (
    BYTES_CODEC,
    NUMPY_CODEC,
    OPENFF_QUANTITY_CODEC,
    OPENFF_UNIT_CODEC,
    PATH_CODEC,
    SETTINGS_CODEC,
)
from gufe.custom_json import JSONSerializerDeserializer, custom_json_factory
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


class CustomJSONCodingTest:
    """Base class for testing codecs.

    In ``setup()``, user must define the following:

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
        self.objs = [np.array([[1.0, 0.0], [2.0, 3.2]]), np.array([1, 0])]
        shapes = [[2, 2], [2,]]
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
            json_str_2 = json.dumps(obj, cls=encoder)
            assert json_str == json_str_2


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
    def setup(self):
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
                    "remove_com": False,
                    "hydrogen_mass": 3.0,
                    "forcefields": [
                        "amber/ff14SB.xml",
                        "amber/tip3p_standard.xml",
                        "amber/tip3p_HFE_multivalent.xml",
                        "amber/phosaa10.xml",
                    ],
                    "small_molecule_forcefield": "openff-2.0.0",
                },
                "thermo_settings": {
                    "__class__": "ThermoSettings",
                    "__module__": "gufe.settings.models",
                    ":is_custom:": True,
                    "temperature": {
                        "magnitude": 300,
                        "unit": "kelvin",
                        ":is_custom:": True,
                        "pint_unit_registry": "openff_units",
                    },
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


class TestOpenFFQuanityCodec(CustomJSONCodingTest):
    def setup(self):
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


class TestOpenFFUnitCodec(CustomJSONCodingTest):
    def setup(self):
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
