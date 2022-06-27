# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import abc
import hashlib
import pathlib
import shutil
import io
import os
import glob
from typing import Union, Tuple, ContextManager
import dataclasses

from ..errors import (
    MissingExternalResourceError, ChangedExternalResourceError
)


@dataclasses.dataclass
class Metadata:
    # as a dataclass to facilitate inheritance and type checking:
    # https://stackoverflow.com/a/50369898
    md5: str

    def to_dict(self):
        return dataclasses.asdict(self)


class _ForceContext:
    """Wrapper that forces objects to only be used a context managers.

    Filelike objects can often be used with explicit open/close. This
    requires the returned byteslike to be consumed as a context manager.
    """
    def __init__(self, context):
        self._context = context

    def __enter__(self):
        return self._context.__enter__()

    def __exit__(self, exc_type, exc_value, traceback):
        return self._context.__exit__(exc_type, exc_value, traceback)


class ExternalStorage(abc.ABC):
    """Abstract base for external storage."""

    def _get_hexdigest(self, location) -> str:
        """Default method for getting the md5 hexdigest.

        Subclasses may implement faster approaches.
        """
        with self._load_stream(location) as byteslike:
            hasher = hashlib.md5()
            # TODO: chunking may give better performance
            hasher.update(byteslike.read())
            digest = hasher.hexdigest()

        return digest

    def get_metadata(self, location: str) -> Metadata:
        """
        Obtain the metadata associated with the actual stored data.

        We always and only obtain the metadata *after* the data has been
        stored. This is because some potential metadata fields, such as
        last-modified timestamps, may not be known until the data is stored.

        Parameters
        ----------
        location : str
            the label to obtain the metadata about

        Returns
        -------
        Metadata :
            Metadata for this object.
        """
        # NOTE: in the future, this may become a (named)tuple of metadata.
        # Subclasses would implement private methods to get each field.
        return Metadata(md5=self._get_hexdigest(location))

    def get_filename(self, location) -> str:
        # we'd like to not need to include the get_filename method, but for
        # now some consumers of this may not work with byteslike
        return self._get_filename(location)

    def load_stream(self, location) -> ContextManager:
        """
        Load data for the given chunk in a context manager.

        This returns a ``_ForceContext``, which requires that the returned
        object be used as a context manager. That ``_ForcedContext`` should
        wrap a byteslike objet.

        Subclasses should implement ``_load_stream``.

        Parameters
        ----------
        location : str
            the label for the data to load
        metadata : str
            metadata to validate that the loaded data is still valid

        Returns
        _ForceContext :
            Wrapper around the byteslike
        """
        return _ForceContext(self._load_stream(location))

    def delete(self, location):
        """
        Delete an existing data chunk from the backend.

        Subclasses should implement ``_delete``.

        Parameters
        ----------
        location : str
            label for the data to delete

        Raises
        ------
        MissingExternalResourceError
            If the resource to be deleted does not exist
        """
        return self._delete(location)

    def store_bytes(self, location, byte_data):
        """
        Store given data in the backend.

        Subclasses should implement ``_store``.

        Parameters
        ----------
        location : str
            label associated with the data to store
        byte_data : bytes
            bytes to store
        """
        return self._store_bytes(location, byte_data)

    def store_path(self, location, path: pathlib.Path):
        self._store_path(location, path)

    def exists(self, location) -> bool:
        """
        Check whether a given label has already been used.

        Subclasses should implement ``_exists``.

        Parameters
        ----------
        location : str
            the label to check for

        Return
        ------
        bool
            True if this label has associated data in the backend, False if
            not
        """
        return self._exists(location)

    def iter_contents(self, prefix=""):
        """Iterate over the labels in this storage.

        Parameters
        ----------
        prefix : str
            Only iterate over paths that start with the given prefix.

        Returns
        -------
        Iterator[str] :
            Contents of this storage, which may include items without
            metadata.
        """
        return self._iter_contents(prefix)

    @abc.abstractmethod
    def _iter_contents(self, prefix=""):
        raise NotImplementedError()

    @abc.abstractmethod
    def _store_bytes(self, location, byte_data):
        """
        For implementers: This should be blocking, even if the storage
        backend allows asynchronous storage.
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def _store_path(self, location, path):
        """
        For implementers: This should be blocking, even if the storage
        backend allows asynchronous storage.
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def _exists(self, location) -> bool:
        raise NotImplementedError()

    @abc.abstractmethod
    def _delete(self, location):
        raise NotImplementedError()

    @abc.abstractmethod
    def _get_filename(self, location):
        raise NotImplementedError()

    @abc.abstractmethod
    def _load_stream(self, location):
        raise NotImplementedError()
