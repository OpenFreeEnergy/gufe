import abc
import hashlib
import warnings
import pathlib
import io
import contextlib
import functools

from typing import ClassVar, Union

from gufe.storage.errors import (
    MissingExternalResourceError, ChangedExternalResourceError
)


class _ForceContext:
    """Wrapper that forces objects to only be used a context managers.

    Filelike objects can often be used with explicit open/close. This
    requires the returned filelike to be consumed as a context manager.
    """
    def __init__(self, context):
        self._context = context

    def __enter__(self):
        return self._context.__enter__()

    def __exit__(self, exc_type, exc_value, traceback):
        self._context.__exit__(exc_type, exc_value, traceback)


class ExternalStorage(abc.ABC):
    """Abstract base for external storage."""

    allow_changed: ClassVar[bool] = False

    def validate(self, location, metadata):
        if not self.get_metadata(location) == metadata:
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


    def get_metadata(self, location: str):
        with self._load_stream(location, bytes_mode=True) as filelike:
            hasher = hashlib.md5()
            # TODO: chunking may give better performance
            hasher.update(filelike.read())
            digest =  hasher.hexdigest()

        return digest


    def get_filename(self, location, metadata) -> str:
        # we'd like to not need to include the get_filename method, but for
        # now some consumers of this may not work with filelike
        self.validate(location, metadata)
        return self._get_filename(location)

    # @force_context
    def load_stream(self, location, metadata, bytes_mode: bool = True):
        self.validate(location, metadata)
        return _ForceContext(self._load_stream(location, bytes_mode))

    @abc.abstractmethod
    def store(self, location, byte_data):
        raise NotImplementedError()

    @abc.abstractmethod
    def exists(self, location):
        raise NotImplementedError()

    @abc.abstractmethod
    def delete(self, location):
        raise NotImplementedError()

    @abc.abstractmethod
    def _get_filename(self, location):
        raise NotImplementedError()

    @abc.abstractmethod
    def _load_stream(self, location, bytes_mode: bool = True):
        raise NotImplementedError()


# TODO: this should use pydantic to check init inputs
class FileStorage(ExternalStorage):
    def __init__(self, root_dir: Union[pathlib.Path, str]):
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

        return str(path), self.get_metadata(path)

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

    def _get_filename(self, location):
        return str(self._as_path(location))

    def _load_stream(self, location, bytes_mode: bool = True):
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
        self._data[location] = byte_data
        return location, self.get_metadata(location)

    def _get_filename(self, location):
        # TODO: how to get this to work? how to manage tempfile? maybe a
        # __del__ here?
        pass

    def _load_stream(self, location, bytes_mode: bool = True):
        byte_data = self._data[location]
        if bytes_mode:
            stream = io.BytesIO(byte_data)
        else:
            # I guess we just have to assume UTF in this case
            stream = io.StringIO(byte_data.decode('utf-8'))

        return stream
