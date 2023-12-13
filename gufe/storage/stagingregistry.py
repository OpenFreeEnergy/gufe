from __future__ import annotations

from typing import Union, Optional
from pathlib import Path
from os import PathLike, rmdir, remove
from .externalresource import ExternalStorage, FileStorage
from contextlib import contextmanager

from gufe.utils import delete_empty_dirs

import logging
_logger = logging.getLogger(__name__)

def _safe_to_delete_file(
    external: ExternalStorage,
    path: PathLike
) -> bool:
    """Check that deleting this file will not remove it from external"""
    # kind of brittle: deals with internals of FileStorage
    if isinstance(external, FileStorage):
        root = external.root_dir
    else:
        return True

    p = Path(path)
    try:
        label = str(p.relative_to(root))
    except ValueError:
        return True
    return not external.exists(label)


class StagingRegistry:
    """PathLike local representation of an :class:`.ExternalStorage`.

    This connects objects on a local filesystem to the key-value store of a
    (possibly remote) :class:`.ExternalStorage`. It presents a FileLike
    interface to users, but internally (via the :class:`.StagingPath` objects
    it contains in its registry) maps local filenames to the keys (labels)
    for the key-value store.

    1. If a local path is requested that corresponds to an existing label in
       the :class:`.ExternalStorage`, this object will "download" the
       contents of that key to that local path.

    2. When requested, it transfers any newly created files to the
       :class:`.ExternalStorage`.

    3. It can delete all of the files it manages.

    Parameters
    ----------
    scratch : PathLike
        the scratch directory shared by all objects on this host
    external : :class:`.ExternalStorage`
        external storage resource where objects should eventualy go
    staging : PathLike
        name of the subdirectory of scratch where staged results are
        temporarily stored; default is '.staging'. This must be the same for
        all units within a DAG.
    delete_staging : bool
        whether to delete the contents of the $SCRATCH/$STAGING
        directory when this object is deleted
    """
    def __init__(
        self,
        scratch: PathLike,
        external: ExternalStorage,
        *,
        staging: PathLike = Path(".staging"),
        delete_staging: bool = True,
        keep_empty_dirs: bool = False,
    ):
        self.external = external
        self.scratch = Path(scratch)
        self.delete_staging = delete_staging
        self.keep_empty_dirs = keep_empty_dirs
        self.staging = staging

        self.registry: set[StagingPath] = set()
        self.preexisting: set[StagingPath] = set()
        self.staging_dir = self.scratch / staging
        self.staging_dir.mkdir(exist_ok=True, parents=True)

    def _delete_file_safe(self, file):
        """Check if deleting this file will remove it from external."""
        return _safe_to_delete_file(
            external=self.external,
            path=file
        )

    def transfer_single_file_to_external(self, held_file: StagingPath):
        """Transfer a given file from staging into external storage
        """
        path = held_file.as_path()
        if not path.exists():
            _logger.info(f"Found nonexistent path {path}, not "
                         "transfering to external storage")
        elif path.is_dir():
            _logger.debug(f"Found directory {path}, not "
                          "transfering to external storage")
        else:
            _logger.info(f"Transfering {path} to external storage")
            self.external.store_path(held_file.label, path)
            return held_file

        return None  # no transfer

    def transfer_staging_to_external(self):
        """Transfer all objects in the registry to external storage

        """
        return [
            transferred
            for file in self.registry
            if (transferred := self.transfer_single_file_to_external(file))
        ]

    def _delete_file(self, file: StagingPath):
        path = file.as_path()
        if path.exists():
            if not path.is_dir():
                _logger.debug(f"Removing file '{file}'")
                path.unlink()
            else:
                _logger.debug(
                    f"During staging cleanup, the directory '{file}' was "
                    "found as a staged path. This will be deleted only if "
                    "`keep_empty` is False."
                )
            self.registry.remove(file)
        else:
            _logger.warning(
                f"During staging cleanup, file '{file}' was marked for "
                "deletion, but can not be found on disk."
            )

    def cleanup(self):
        """Perform end-of-lifecycle cleanup.
        """
        if self.delete_staging:
            for file in self.registry - self.preexisting:
                if self._delete_file_safe(file):
                    self._delete_file(file)

            if not self.keep_empty_dirs:
                delete_empty_dirs(self.staging_dir)

    def register_path(self, staging_path: StagingPath):
        """
        Register a :class:`.StagingPath` with this :class:`.StagingRegistry`.

        This marks a given path as something for this object to manage, by
        loading it into the ``registry``. This way it is tracked such that
        its contents can be transfered to the :class:`.ExternalStorage` and
        such that the local copy can be deleted when it is no longer needed.

        If this objects's :class:`.ExternalStorage` already has data for the
        label associated with the provided :class:`.Stagingpath`, then the
        contents of that should copied to the local path so that it can be
        read by the user.

        Parameters
        ----------
        staging_path: :class:`.StagingPath`
            the path to track
        """
        label_exists = self.external.exists(staging_path.label)
        fspath = staging_path.as_path()

        # TODO: what if the staging path is a directory? not sure that we
        # have a way to know that; but not sure that adding it to the
        # registry is right either
        if not fspath.parent.exists():
            fspath.parent.mkdir(parents=True, exist_ok=True)

        self.registry.add(staging_path)

        # if this is a file that exists, bring it into our subdir
        # NB: this happens even if you're intending to overwrite the path,
        # which is kind of wasteful
        if label_exists:
            self._load_file_from_external(self.external, staging_path)

    def _load_file_from_external(self, external: ExternalStorage,
                                 staging_path: StagingPath):
        # import pdb; pdb.set_trace()
        scratch_path = self.staging_dir / staging_path.path
        # TODO: switch this to using `get_filename` and `store_path`
        if scratch_path.exists():
            self.preexisting.add(staging_path)

        with external.load_stream(staging_path.label) as f:
            external_bytes = f.read()
            ...  # TODO: check that the bytes are the same if preexisting?

        scratch_path.parent.mkdir(exist_ok=True, parents=True)
        with open(scratch_path, mode='wb') as f:
            f.write(external_bytes)

    def __truediv__(self, path: Union[PathLike, str]):
        return StagingPath(root=self, path=path)

    def __fspath__(self):
        return str(self.staging_dir)

    def __repr__(self):
        return (
            f"{self.__class__.__name__}('{self.scratch}', {self.external})"
        )

    def __del__(self):  # -no-cov-
        # in case someone doesn't use this within a context manager
        if self.staging_dir.exists():
            self.cleanup()


class SharedStaging(StagingRegistry):
    """Staging for shared external storage.

    This enables read-only versions to be loaded from other units.
    """
    def __init__(
        self,
        scratch: PathLike,
        external: ExternalStorage,
        *,
        staging: PathLike = Path(".staging"),
        delete_staging: bool = True,
        keep_empty_dirs: bool = False,
        read_only: bool = False,
    ):
        super().__init__(scratch, external, staging=staging,
                         delete_staging=delete_staging,
                         keep_empty_dirs=keep_empty_dirs)
        self.read_only = read_only

    def transfer_single_file_to_external(self, held_file: StagingPath):
        if self.read_only:
            _logger.debug("Read-only: Not transfering to external storage")
            return  # early exit

        return super().transfer_single_file_to_external(held_file)

    def transfer_staging_to_external(self):
        if self.read_only:
            _logger.debug("Read-only: Not transfering to external storage")
            return  # early exit

        return super().transfer_staging_to_external()

    def register_path(self, staging_path: StagingPath):
        label_exists = self.external.exists(staging_path.label)

        if self.read_only and not label_exists:
            raise IOError(f"Unable to create '{staging_path.label}'. File "
                          "does not exist in external storage, and this "
                          "staging path is read-only.")

        super().register_path(staging_path)


class PermanentStaging(StagingRegistry):
    """Staging directory for the permanent storage.

    This allows files to be downloaded from a shared
    :class:`.ExternalStorage`.
    """
    def __init__(
        self,
        scratch: PathLike,
        external: ExternalStorage,
        shared: ExternalStorage,
        *,
        staging: PathLike = Path(".staging"),
        delete_staging: bool = True,
        keep_empty_dirs: bool = False,
    ):
        super().__init__(scratch, external, staging=staging,
                         delete_staging=delete_staging,
                         keep_empty_dirs=keep_empty_dirs)
        self.shared = shared

    def _delete_file_safe(self, file):
        shared_safe = _safe_to_delete_file(
            external=self.shared,
            path=file
        )
        return shared_safe and super()._delete_file_safe(file)

    def transfer_single_file_to_external(self, held_file: StagingPath):
        # if we can't find it locally, we load it from shared storage
        path = held_file.as_path()
        if not path.exists():
            self._load_file_from_external(self.shared, held_file)

        super().transfer_single_file_to_external(held_file)


class StagingPath:
    """PathLike object linking local path with label for external storage.

    On creation, this registers with a :class:`.StagingRegistry` that will
    manage the local path and transferring data with its
    :class:`.ExternalStorage`.

    This object can always be used as a FileLike (using, e.g., the standard
    ``open`` builtin). This requires that a staged path that exists on an
    external resource be downloaded into a local file when it is referenced.

    For a representation of a file that does not require the download (for
    example, when deserializing results that point to files) instead use
    :class:`.ExternalFile`.
    """
    def __init__(self, root: StagingRegistry,
                 path: Union[PathLike, str]):
        self.root = root
        self.path = Path(path)

    def register(self):
        """Register this path with its StagingRegistry.

        If a file associated with this path exists in an external storage,
        it will be downloaded to the staging area as part of registration.
        """
        self.root.register_path(self)

    def __truediv__(self, path: Union[PathLike, str]):
        return StagingPath(self.root, self.path / path)

    def __eq__(self, other):
        return (isinstance(other, StagingPath)
                and self.root == other.root
                and self.path == other.path)

    def __hash__(self):
        return hash((self.root, self.path))

    def as_path(self):
        """Return the pathlib.Path where this is staged"""
        return Path(self._fspath)

    @property
    def _fspath(self):
        return str(self.root.staging_dir / self.path)

    def __fspath__(self):
        self.register()
        return self._fspath

    @property
    def label(self) -> str:
        """Label used in :class:`.ExternalStorage` for this path"""
        return str(self.path)

    def __repr__(self):
        return f"StagingPath('{self._fspath}')"

    # TODO: how much of the pathlib.Path interface do we want to wrap?
    # although edge cases may be a pain, we can get most of it with, e.g.:
    # def exists(self): return Path(self).exists()
    # but also, can do pathlib.Path(staging_path) and get hte whole thing
