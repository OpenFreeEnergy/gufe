import pytest
from gufe.molhashing import serialize_numpy, deserialize_numpy
import numpy as np


@pytest.mark.parametrize('dtype', ['float32', 'float64', 'int'])
def test_numpy_serialization_cycle(dtype):
    arr = np.array([1, 2], dtype=dtype)
    ser = serialize_numpy(arr)
    deser = deserialize_numpy(ser)
    reser = serialize_numpy(deser)

    np.testing.assert_equal(arr, deser)
    assert ser == reser


def test_numpy_serialize_different_dtypes():
    arr32 = np.array([1, 2], dtype=np.float32)
    arr64 = np.array([1, 2], dtype=np.float64)
    assert serialize_numpy(arr32) != serialize_numpy(arr64)
