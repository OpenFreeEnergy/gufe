import hashlib
import warnings
import pathlib

from typing import ClassVar

from gufe.storage.errors import (
    MissingExternalResourceError, ChangedExternalResourceError
)


class ExternalStorage:
    """Abstract base for external storage."""

    allow_changed: ClassVar[bool] = False

    def validate(self, location, sha2):
        if not self.get_sha2(location) == sha2:
            msg = (f"Hash mismatch for {location}: this object "
                   "may have changed.")
            if not self.allow_changed:
                # NOTE: having it here instead of in a DataLoader means that
                # you can only change it for the whole system, instead of
                # for a specific object
                raise ChangedExternalResourceError(
                    msg + " To allow this, set ExternalStorage."
                    "allow_changed = True"
                )
            else:
                warnings.warn(msg)

    def get_sha2(self, location: str):
        filelike = self._as_filelike(location, bytes_mode=True)
        sha2 = hashlib.sha256()
        # TODO: this can be probably be improved by chunking
        sha2.update(filelike.read())
        return sha2.hexdigest()

    def as_filename(self, location, sha2) -> str:
        # we'd like to not need to include the as_filename method, but for
        # now some consumers of this may not work with filelike
        self.validate(location, sha2)
        return self._as_filename(location)

    def as_filelike(self, location, sha2, bytes_mode: bool = True):
        self.validate(location, sha2)
        return self._as_filelike(location, bytes_mode)

    def store(self, location, byte_data):
        raise NotImplementedError()

    def exists(self, location):
        raise NotImplementedError()

    def delete(self, location):
        raise NotImplementedError()

    def _as_filename(self, location):
        raise NotImplementedError()

    def _as_filelike(self, location, bytes_mode: bool = True):
        raise NotImplementedError()


# TODO: this should use pydantic to check init inputs
class FileStorage(ExternalStorage):
    def __init__(self, root_dir):
        self.root_dir = pathlib.Path(root_dir)

    def exists(self, location):
        return self._as_path(location).exists()

    def store(self, location, byte_data):
        path = self._as_path(location)
        directory = path.parent
        filename = path.name
        directory.mkdir(parents=True, exist_ok=True)
        with open(path, mode='wb') as f:
            f.write(byte_data)

        return (str(path), self.get_sha2(path))

    def delete(self, location):
        path = self._as_path(location)
        if self.exists(location):
            path.unlink()
        else:
            raise MissingExternalResourceError(
                f"Unable to delete f{str(path)}: File does not exist"
            )

    def _as_path(self, location):
        return self.root_dir / pathlib.Path(location)

    def _as_filename(self, location):
        return str(self._as_path(location))

    def _as_filelike(self, location, bytes_mode: bool = True):
        mode = 'r'
        if bytes_mode:
            mode += 'b'

        try:
            return open(self._as_path(location), mode)
        except OSError as e:
            raise MissingExternalResourceError(str(e))


class MemoryStorage(ExternalStorage):
    """Not for production use, but potentially useful in testing"""
    def __init__(self):
        self._data = {}

    def exists(self, location):
        return location in self._data

    def store(self, location, byte_data):
        sha2 = hashlib.sha256()
        sha2.update(byte_data)
        self._data[location] = byte_data
        return location, sha2.hexdigest()

    def _as_filename(self, location):
        # TODO: how to get this to work? how to manage tempfile? maybe a
        # __del__ here?
        pass

    def _as_filelike(self, location, bytes_mode: bool = True):
        byte_data = self._data[location]
        if bytes_mode:
            stream = io.BytesIO(byte_data)
        else:
            # I guess we just have to assume UTF in this case
            stream = io.StringIO(byte_data.encode('utf-8'))

        return stream
