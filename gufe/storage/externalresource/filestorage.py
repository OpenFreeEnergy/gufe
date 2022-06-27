# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import pathlib
import shutil
import os
from typing import Union, Tuple, ContextManager

from .base import ExternalStorage

from ..errors import (
    MissingExternalResourceError, ChangedExternalResourceError
)


# TODO: this should use pydantic to check init inputs
class FileStorage(ExternalStorage):
    def __init__(self, root_dir: Union[pathlib.Path, str]):
        self.root_dir = pathlib.Path(root_dir)

    def _exists(self, location):
        return self._as_path(location).exists()

    def _store_bytes(self, location, byte_data):
        path = self._as_path(location)
        directory = path.parent
        filename = path.name
        # TODO: add some stuff here to catch permissions-based errors
        directory.mkdir(parents=True, exist_ok=True)
        with open(path, mode='wb') as f:
            f.write(byte_data)

    def _store_path(self, location, path):
        my_path = self._as_path(location)
        if path.resolve() != my_path.resolve():
            shutil.copyfile(path, my_path)

    def _iter_contents(self, prefix):
        start_dir = (self.root_dir / pathlib.Path(prefix).parent).resolve()
        for dirpath, _, filenames in os.walk(start_dir):
            for filename in filenames:
                path = pathlib.Path(dirpath) / filename
                location = self._get_location(path)
                if location.startswith(prefix):
                    yield location

    def _delete(self, location):
        path = self._as_path(location)
        if self.exists(location):
            path.unlink()
        else:
            raise MissingExternalResourceError(
                f"Unable to delete '{str(path)}': File does not exist"
            )

    def _as_path(self, location):
        return self.root_dir / pathlib.Path(location)

    def _get_location(self, path: Union[str, pathlib.Path]) -> str:
        # this is essentially the reverse behavior of _as_path
        path = pathlib.Path(path)
        relpath = path.relative_to(self.root_dir.resolve())
        return str(relpath)

    def _get_filename(self, location):
        return str(self._as_path(location))

    def _load_stream(self, location):
        try:
            return open(self._as_path(location), 'rb')
        except OSError as e:
            raise MissingExternalResourceError(str(e))
