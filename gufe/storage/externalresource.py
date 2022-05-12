import hashlib
import warnings
import pathlib

from typing import ClassVar


class ExternalResourceError(Exception):
    """Base class for errors due to problems with external resources"""
    # TODO: is it necessary to have a base class here? Would you ever have
    # one catch that handles both subclass errors?


class MissingExternalResourceError(ExternalResourceError):
    """Error when the external resource could not be loaded"""


class ChangedExternalResourceError(ExternalResourceError):
    """Error when there's a SHA256 mismatch with an external resource"""


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
        return self._as_filename(self.location)

    def as_filelike(self, location, sha2, bytes_mode: bool = True):
        self.validate(location, sha2)
        return self._as_filelike(location, bytes_mode)

    def create(self, location, byte_data):
        raise NotImplementedError()

    def exists(self, location):
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

    def create(self, location, byte_data):
        path = self._as_path(location)
        directory = path.directory
        filename = path.name
        os.makedir(directory, existsok=True)
        with open(path, mode='wb') as f:
            f.write(byte_data)

        return (str(path), cls.get_sha2(path))

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
