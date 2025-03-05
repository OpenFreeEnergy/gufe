import abc
import uuid
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import openff
import pytest
from numpy import testing as npt

from gufe.serialization.msgpack import packb, unpackb
from gufe.settings import models


def test_unpack_unknown_extension():
    payload = bytes([0xD4, 0x7F, 0x00])
    with pytest.raises(ValueError, match="Found an unknown extension code: 127"):
        unpackb(payload)


class CustomMessagePackCodingTest(abc.ABC):

    @abc.abstractmethod
    def setup_method(self):
        return NotImplementedError

    def _equality_check(self, original, reconstructed):
        assert original == reconstructed

    def test_round_trip(self):
        for obj in self.objs:
            msgpack_bytes = packb(obj)
            reconstructed = unpackb(msgpack_bytes)

            self._equality_check(obj, reconstructed)


class TestNumpyCoding(CustomMessagePackCodingTest):

    def setup_method(self):
        self.objs = [
            np.array([[1.0, 0.0], [2.0, 3.2]]),
            np.array([1, 0]),
            np.array([1.0, 2.0, 3.0], dtype=np.float32),
            np.array([1.0, 2.0, 3.0], dtype=np.float16),
        ]

    def _equality_check(self, original, reconstructed):
        npt.assert_array_equal(original, reconstructed)
        assert original.dtype == reconstructed.dtype


class TestLargeIntCoding(CustomMessagePackCodingTest):

    def setup_method(self):
        self.objs = [
            -1,
            0,
            1,
            2 ^ 64,
            2 ^ 128,
            2 ^ 256,
        ]


class TestPathCoding(CustomMessagePackCodingTest):

    def setup_method(self):
        self.objs = [
            Path("/"),
            Path("/foo"),
            Path("/foo/bar"),
            Path("."),
            Path("./foo/"),
        ]


class TestNumpyGenericCoding(CustomMessagePackCodingTest):

    def setup_method(self):
        self.objs = [
            np.bool_(True),
            np.float16(1.0),
            np.float32(1.0),
            np.complex128(1.0),
            np.clongdouble(1.0),
            np.uint64(1),
        ]


class TestUnitCoding(CustomMessagePackCodingTest):

    def setup_method(self):
        self.objs = [
            openff.units.DEFAULT_UNIT_REGISTRY("1.0 * kg meter per second squared"),
            1 * openff.units.unit.nanometer,
            1.0 * openff.units.unit.nanometer,
            0.8 * openff.units.unit.nanometer,
            np.array([-1.0, 0, 1.0]) * openff.units.unit.nanometer,
        ]

    def _equality_check(self, original, reconstructed):

        match original.m:
            case np.ndarray():
                npt.assert_array_equal(original, reconstructed)
                assert original.dtype == reconstructed.dtype
            case _:
                assert original == reconstructed


class TestSettingsCoding(CustomMessagePackCodingTest):

    def setup_method(self):
        self.objs = [
            models.Settings.get_defaults(),
        ]


class TestTimeStampCoding(CustomMessagePackCodingTest):

    def setup_method(self):
        self.objs = [datetime.now(timezone.utc)]


class TestUUIDCoding(CustomMessagePackCodingTest):

    def setup_method(self):
        self.objs = [uuid.uuid4() for _ in range(3)]
