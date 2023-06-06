from typing import Union, Optional
from pathlib import Path
from os import PathLike, rmdir, remove
from .externalresource import ExternalStorage, FileStorage
from contextlib import contextmanager

import logging
_logger = logging.getLogger(__name__)

def _safe_to_delete_holding(external, path, prefix):
    """Check if deleting ``path`` could delete externally stored data.

    If external storage is a FileStorage, then it will storage files for
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


def _delete_empty_dirs(root, delete_root=True):
    """Delete all empty directories.

    Repeats so that directories that only contained empty directories also
    get deleted.
    """
    root = Path(root)

    def find_empty_dirs(directory):
        if not (paths := list(directory.iterdir())):
            return [directory]
        directories = [p for p in paths if p.is_dir()]
        return sum([find_empty_dirs(d) for d in directories], [])


    while root.exists() and (empties := find_empty_dirs(root)):
        if empties == [root] and not delete_root:
            return
        for directory in empties:
            _logger.debug(f"Removing '{directory}'")
            rmdir(directory)


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

    3. It can delete all of the files it manages

    This can be opened in "read-only" mode, which prevents new files from
    being created, but does not prevent changes to existing versions of
    local files.

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
    holding : PathLike
        name of the subdirectory of scratch where staged results are
        temporarily stored; default is '.holding'. This must be the same for
        all units within a DAG.
    delete_holding : bool
        whether to delete the contents of the $SCRATCH/$HOLDING/$PREFIX
        directory when this object is deleted
    read_only : bool
        write to prevent NEW files from being written within this staging
        directory. NOTE: This will not prevent overwrite of existing files
        in scratch space, but it will prevent changed files from uploading
        to the external storage.
    """
    def __init__(
        self,
        scratch: PathLike,
        external: ExternalStorage,
        prefix: str,
        *,
        holding: PathLike = Path(".holding"),
        delete_holding: bool = True,
    ):
        self.external = external
        self.scratch = Path(scratch)
        self.prefix = Path(prefix)
        self.delete_holding = delete_holding
        self.holding = holding

        self.registry : set[StagingPath] = set()
        self.preexisting : set[StagingPath] = set()
        self.staging_dir = self.scratch / holding / prefix
        self.staging_dir.mkdir(exist_ok=True, parents=True)

    def _delete_holding_safe(self):
        return _safe_to_delete_holding(
            external=self.external,
            path=self.staging_dir,
            prefix=self.prefix,
        )

    def transfer_single_file_to_external(self, held_file):
        """Transfer a given file from holding into external storage
        """

        path = Path(held_file)
        if not path.exists():
            logging.info(f"Found nonexistent path {path}, not "
                         "transfering to external storage")
        elif path.is_dir():
            logging.debug(f"Found directory {path}, not "
                          "transfering to external storage")
        else:
            logging.info(f"Transfering {path} to external storage")
            self.external.store_path(held_file.label, path)

    def transfer_holding_to_external(self):
        """Transfer all objects in the registry to external storage"""
        for obj in self.registry:
            self.transfer_single_file_to_external(obj)

    def cleanup(self):
        """Perform end-of-lifecycle cleanup.
        """
        if self.delete_holding and self._delete_holding_safe():
            for file in self.registry - self.preexisting:
                remove(file)
            _delete_empty_dirs(self.staging_dir)

    def register_path(self, staging_path):
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

        self.registry.add(staging_path)

        # if this is a file that exists, bring it into our subdir
        # NB: this happens even if you're intending to overwrite the path,
        # which is kind of wasteful
        if label_exists:
            self._load_file_from_external(self.external, staging_path)

    def _load_file_from_external(self, external, staging_path):
            scratch_path = self.staging_dir / staging_path.path
            # TODO: switch this to using `get_filename` and `store_path`
            with external.load_stream(staging_path.label) as f:
                external_bytes = f.read()
            if scratch_path.exists():
                self.preexisting.add(staging_path)
                ... # TODO: something to check that the bytes are the same?
            scratch_path.parent.mkdir(exist_ok=True, parents=True)
            with open(scratch_path, mode='wb') as f:
                f.write(external_bytes)

    def __truediv__(self, path: PathLike):
        return StagingPath(root=self, path=path)

    def __fspath__(self):
        return str(self.staging_dir)

    def __repr__(self):
        return (
            f"StagingDirectory({self.scratch}, {self.external}, "
            f"{self.prefix})"
        )


class SharedStaging(StagingDirectory):
    def __init__(
        self,
        scratch: PathLike,
        external: ExternalStorage,
        prefix: str,
        *,
        holding: PathLike = Path(".holding"),
        delete_holding: bool = True,
        read_only: bool = False,
    ):
        super().__init__(scratch, external, prefix, holding=holding,
                         delete_holding=delete_holding)
        self.read_only = read_only

    def get_other_shared(self, prefix, delete_holding=None):
        """Get a related unit's staging directory.
        """
        if delete_holding is None:
            delete_holding = self.delete_holding

        return SharedStaging(
            scratch=self.scratch,
            external=self.external,
            prefix=prefix,
            holding=self.holding,
            delete_holding=delete_holding,
            read_only=True,
        )

    @contextmanager
    def other_shared(self, prefix, delete_holding=None):
        """Context manager approach for getting a related unit's directory.

        This is usually the recommended way to get a previous unit's shared
        data.
        """
        other = self.get_other_shared(prefix, delete_holding)
        yield other
        other.cleanup()

    def transfer_single_file_to_external(self, held_file):
        if self.read_only:
            logging.debug("Read-only: Not transfering to external storage")
            return  # early exit

        super().transfer_single_file_to_external(held_file)

    def transfer_holding_to_external(self):
        if self.read_only:
            logging.debug("Read-only: Not transfering to external storage")
            return  # early exit

        super().transfer_holding_to_external()

    def register_path(self, staging_path):
        label_exists = self.external.exists(staging_path.label)

        if self.read_only and not label_exists:
            raise IOError(f"Unable to create '{staging_path.label}'. File "
                          "does not exist in external storage, and This "
                          "staging path is read-only.")

        super().register_path(staging_path)


class PermanentStaging(StagingDirectory):
    def __init__(
        self,
        scratch: PathLike,
        external: ExternalStorage,
        shared: ExternalStorage,
        prefix: str,
        *,
        holding: PathLike = Path(".holding"),
        delete_holding: bool = True,
    ):
        super().__init__(scratch, external, prefix, holding=holding,
                         delete_holding=delete_holding)
        self.shared = shared

    def _delete_holding_safe(self):
        shared_safe = _safe_to_delete_holding(
            external=self.shared,
            path=self.staging_dir,
            prefix=self.prefix
        )
        return shared_safe and super()._delete_holding_safe()

    def transfer_single_file_to_external(self, held_file):
        # if we can't find it locally, we load it from shared storage
        path = Path(held_file)
        if not path.exists():
            self._load_file_from_external(self.shared, held_file)

        super().transfer_single_file_to_external(held_file)


class StagingPath:
    """PathLike object linking local path with label for external storage.

    On creation, this registers with a :class:`.StagingDirectory` that will
    manage the local path and transferring data with its
    :class:`.ExternalStorage`.
    """
    def __init__(self, root: StagingDirectory, path: PathLike):
        self.root = root
        self.path = Path(path)
        self.root.register_path(self)

    def __truediv__(self, path):
        return StagingPath(self.root, self.path / path)

    def __fspath__(self):
        return str(self.root.staging_dir / self.path)

    @property
    def label(self):
        """Label used in :class:`.ExternalStorage` for this path"""
        return str(self.root.prefix / self.path)

    def __repr__(self):
        return f"StagingPath({self.__fspath__()})"

    # TODO: how much of the pathlib.Path interface do we want to wrap?
    # although edge cases may be a pain, we can get most of it with, e.g.:
    # def exists(self): return Path(self).exists()
    # but also, can do pathlib.Path(staging_path) and get hte whole thing
