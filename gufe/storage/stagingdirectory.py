from __future__ import annotations

from typing import Union, Optional
from pathlib import Path
from os import PathLike, rmdir, remove
from .externalresource import ExternalStorage, FileStorage
from contextlib import contextmanager

from gufe.utils import delete_empty_dirs

import logging
_logger = logging.getLogger(__name__)


def _safe_to_delete_staging(external: ExternalStorage, path: PathLike,
                            prefix: Union[PathLike, str]) -> bool:
    """Check if deleting ``path`` could delete externally stored data.

    If external storage is a FileStorage, then it will store files for
    this unit or dag in the directory ``external.root_dir / prefix``, where
    ``prefix`` is either the unit label or the dag label. If ``path`` is
    inside that directory, then deleting it may delete information from the
    external storage. In that case, this returns False, indicating a
    conflict. Otherwise, this returns True.
    """
    # this is a little brittle; I don't like hard-coding the class here
    if isinstance(external, FileStorage):
        root = Path(external.root_dir) / prefix
    else:
        return True

    p = Path(path)
    try:
        _ = p.relative_to(root)
    except ValueError:
        return True
    else:
        return False



class StagingDirectory:
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
    prefix : str
        label for this specific unit; this should be a slash-separated
        description of where this unit fits in the hierarchy. For example,
        it might be ``$DAG_LABEL/$UNIT_LABEL`` or
        ``$DAG_LABEL/$UNIT_LABEL/$UNIT_REPEAT``. It must be a unique
        identifier for this unit within the permanent storage.
    staging : PathLike
        name of the subdirectory of scratch where staged results are
        temporarily stored; default is '.staging'. This must be the same for
        all units within a DAG.
    delete_staging : bool
        whether to delete the contents of the $SCRATCH/$HOLDING/$PREFIX
        directory when this object is deleted
    """
    def __init__(
        self,
        scratch: PathLike,
        external: ExternalStorage,
        prefix: str,
        *,
        staging: PathLike = Path(".staging"),
        delete_staging: bool = True,
        delete_empty_dirs: bool = True,
    ):
        self.external = external
        self.scratch = Path(scratch)
        self.prefix = Path(prefix)
        self.delete_staging = delete_staging
        self.delete_empty_dirs = delete_empty_dirs
        self.staging = staging

        self.registry : set[StagingPath] = set()
        self.preexisting : set[StagingPath] = set()
        self.staging_dir = self.scratch / staging / prefix
        self.staging_dir.mkdir(exist_ok=True, parents=True)

    def _delete_staging_safe(self):
        """Check if deleting staging will remove data from external.
        """
        return _safe_to_delete_staging(
            external=self.external,
            path=self.staging_dir,
            prefix=self.prefix,
        )

    def transfer_single_file_to_external(self, held_file: StagingPath):
        """Transfer a given file from staging into external storage
        """
        path = Path(held_file.fspath)
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

    def cleanup(self):
        """Perform end-of-lifecycle cleanup.
        """
        if self.delete_staging and self._delete_staging_safe():
            for file in self.registry - self.preexisting:
                path = Path(file.fspath)
                if path.exists():
                    _logger.debug(f"Removing file {file}")
                    path.unlink()
                    self.registry.remove(file)
                else:
                    _logger.warning("During staging cleanup, file "
                                    f"{file} was marked for deletion, but "
                                    "can not be found on disk.")

            if self.delete_empty_dirs:
                delete_empty_dirs(self.staging_dir)

    def register_path(self, staging_path: StagingPath):
        """
        Register a :class:`.StagingPath` with this :class:`.StagingDirectory`.

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
        fspath = Path(staging_path.fspath)
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
            ... # TODO: check that the bytes are the same if preexisting?

        scratch_path.parent.mkdir(exist_ok=True, parents=True)
        with open(scratch_path, mode='wb') as f:
            f.write(external_bytes)

    def __truediv__(self, path: Union[PathLike, str]):
        return StagingPath(root=self, path=path)

    def __fspath__(self):
        return str(self.staging_dir)

    def __repr__(self):
        return (
            f"{self.__class__.__name__}('{self.scratch}', {self.external}, "
            f"{self.prefix})"
        )

    def __del__(self):  # -no-cov-
        # in case someone doesn't use this within a context manager
        if self.staging_dir.exists():
            self.cleanup()


class SharedStaging(StagingDirectory):
    """Staging for shared external storage.

    This enables read-only versions to be loaded from other units.
    """
    def __init__(
        self,
        scratch: PathLike,
        external: ExternalStorage,
        prefix: str,
        *,
        staging: PathLike = Path(".staging"),
        delete_staging: bool = True,
        delete_empty_dirs: bool = True,
        read_only: bool = False,
    ):
        super().__init__(scratch, external, prefix, staging=staging,
                         delete_staging=delete_staging,
                         delete_empty_dirs=delete_empty_dirs)
        self.read_only = read_only

    def _get_other_shared(self, prefix: str,
                         delete_staging: Optional[bool] = None):
        """Get a related unit's staging directory.
        """
        if delete_staging is None:
            delete_staging = self.delete_staging

        return SharedStaging(
            scratch=self.scratch,
            external=self.external,
            prefix=prefix,
            staging=self.staging,
            delete_staging=delete_staging,
            read_only=True,
        )

    @contextmanager
    def other_shared(self, prefix: str,
                     delete_staging: Optional[bool] = None):
        """Context manager approach for getting a related unit's directory.

        This is usually the recommended way to get a previous unit's shared
        data.
        """
        other = self._get_other_shared(prefix, delete_staging)
        yield other
        other.cleanup()

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


class PermanentStaging(StagingDirectory):
    """Staging directory for the permanent storage.

    This allows files to be downloaded from a shared
    :class:`.ExternalStorage`.
    """
    def __init__(
        self,
        scratch: PathLike,
        external: ExternalStorage,
        shared: ExternalStorage,
        prefix: str,
        *,
        staging: PathLike = Path(".staging"),
        delete_staging: bool = True,
        delete_empty_dirs: bool = True,
    ):
        super().__init__(scratch, external, prefix, staging=staging,
                         delete_staging=delete_staging,
                         delete_empty_dirs=delete_empty_dirs)
        self.shared = shared

    def _delete_staging_safe(self):
        shared_safe = _safe_to_delete_staging(
            external=self.shared,
            path=self.staging_dir,
            prefix=self.prefix
        )
        return shared_safe and super()._delete_staging_safe()

    def transfer_single_file_to_external(self, held_file: StagingPath):
        # if we can't find it locally, we load it from shared storage
        path = Path(held_file.fspath)
        if not path.exists():
            self._load_file_from_external(self.shared, held_file)

        super().transfer_single_file_to_external(held_file)


class StagingPath:
    """PathLike object linking local path with label for external storage.

    On creation, this registers with a :class:`.StagingDirectory` that will
    manage the local path and transferring data with its
    :class:`.ExternalStorage`.

    This object can always be used as a FileLike (using, e.g., the standard
    ``open`` builtin). This requires that a staged path that exists on an
    external resource be downloaded into a local file when it is referenced.

    For a representation of a file that does not require the download (for
    example, when deserializing results that point to files) instead use
    :class:`.ExternalFile`.
    """
    def __init__(self, root: StagingDirectory,
                 path: Union[PathLike, str]):
        self.root = root
        self.path = Path(path)

    def __truediv__(self, path: Union[PathLike, str]):
        return StagingPath(self.root, self.path / path)

    @property
    def fspath(self):
        return str(self.root.staging_dir / self.path)

    def __fspath__(self):
        self.root.register_path(self)
        return self.fspath

    @property
    def label(self) -> str:
        """Label used in :class:`.ExternalStorage` for this path"""
        return str(self.root.prefix / self.path)

    def __repr__(self):
        return f"StagingPath('{self.fspath}')"

    # TODO: how much of the pathlib.Path interface do we want to wrap?
    # although edge cases may be a pain, we can get most of it with, e.g.:
    # def exists(self): return Path(self).exists()
    # but also, can do pathlib.Path(staging_path) and get hte whole thing


