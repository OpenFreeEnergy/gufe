# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import io
from typing import Union, Tuple, ContextManager

from .base import ExternalStorage

from ..errors import (
    MissingExternalResourceError, ChangedExternalResourceError
)


class MemoryStorage(ExternalStorage):
    """Not for production use, but potentially useful in testing"""
    def __init__(self):
        self._data = {}

    def _exists(self, location):
        return location in self._data

    def _delete(self, location):
        try:
            del self._data[location]
        except KeyError:
            raise MissingExternalResourceError(
                f"Unable to delete '{location}': key does not exist"
            )

    def _store_bytes(self, location, byte_data):
        self._data[location] = byte_data
        return location, self.get_metadata(location)

    def _store_path(self, location, path):
        with open(path, 'rb') as f:
            byte_data = f.read()

        return self._store_bytes(location, byte_data)

    def _iter_contents(self, prefix):
        for label in self._data:
            if label.startswith(prefix):
                yield label

    def _get_filename(self, location):
        # TODO: how to get this to work? how to manage tempfile? maybe a
        # __del__ here?
        pass

    def _load_stream(self, location):
        byte_data = self._data[location]
        stream = io.BytesIO(byte_data)
        return stream
