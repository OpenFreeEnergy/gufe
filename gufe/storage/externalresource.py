import abc
import hashlib
import pathlib
import shutil
import io
import os
import glob

from typing import Union, Tuple, ContextManager

import boto3

from gufe.storage.errors import (
    MissingExternalResourceError, ChangedExternalResourceError
)


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

    def get_metadata(self, location: str) -> str:
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
        str :
            hexdigest of the md5 hash of the data
        """
        # NOTE: in the future, this may become a (named)tuple of metadata.
        # Subclasses would implement private methods to get each field.
        return self._get_hexdigest(location)

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

    @abc.abstractmethod
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

    def iter_contents(self, prefix):
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

    def iter_contents(self, prefix):
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


class S3Storage(ExternalStorage):
    """File storage backend for AWS S3.

    """
    def __init__(self, session: "boto3.Sessin", bucket: str, prefix: str):
        """

        """
        self.session = session

        self.resource = self.session.resource('s3')
        self.bucket = s3.Bucket(bucket)

        # we don't want to keep any leading or trailing delimiters
        self.prefix = prefix.strip('/')

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
        raise NotImplementedError()

    def _store_bytes(self, location, byte_data):
        """
        For implementers: This should be blocking, even if the storage
        backend allows asynchronous storage.
        """
        key = "/".join([self.prefix, location])

        self.bucket.put_object(Key=key, Body=byte_data)

    def _store_path(self, location, path):
        """
        For implementers: This should be blocking, even if the storage
        backend allows asynchronous storage.
        """
        """
        For implementers: This should be blocking, even if the storage
        backend allows asynchronous storage.
        """
        key = "/".join([self.prefix, location])

        with open(path, 'rb') as f:
            self.bucket.upload_fileobj(f, key)

    def _exists(self, location) -> bool:
        from botocore.exceptions import ClientError

        key = "/".join([self.prefix, location])

        # we do a metadata load as our existence check
        # appears to be most recommended approach
        try:
            self.bucket.Object(key).load()
            return True
        except ClientError:
            return False

    def _delete(self, location):
        key = "/".join([self.prefix, location])

        if self._exists(location):
            self.bucket.Object(key).delete()
        else:
            raise MissingExternalResourceError(
                f"Unable to delete '{str(key)}': Object does not exist"
            )

    def _get_filename(self, location):
        key = "/".join([self.prefix, location])

        object = self.bucket.Object(key)

        url = object.meta.client.generate_presigned_url(
                'get_object',
                ExpiresIn=0,
                Params={'Bucket': self.bucket.name, 'Key': object.key})

        # drop query params from url
        url = url.split('?')[0]

        return url

    def _load_stream(self, location):
        key = "/".join([self.prefix, location])

        try:
            return self.bucket.Object(key).get()['Body']
        except self.resource.meta.client.exceptions.NoSuchKey as e:
            raise MissingExternalResourceError(str(e))
