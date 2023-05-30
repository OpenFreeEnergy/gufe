from typing import Union, Optional
from pathlib import Path
from os import PathLike
from .externalresource import ExternalStorage

import logging
_logger = logging.getLogger(__name__)


def _delete_empty_dirs(root, delete_root=True):
    """Delete all empty directories.

    Repeats so that directories that only contained empty directories also
    get deleted.
    """
    root = Path(root)

    def find_empty_dirs(directory):
        if not (paths := directory.iterdir()):
            return [directory]
        directories = [p for p in paths if p.is_dir()]
        return sum([find_empty_dirs(d) for d in directories], [])

    while empties := find_empty_dirs(root):
        if empties == [root] and not delete_root:
            return
        for directory in empties:
            os.rmdir(directory)


class SharedRoot:
    """PathLike local representation of an :class:`.ExternalStorage`.

    This connects objects on a local filesystem to the key-value store of a
    (possibly remote) :class:`.ExternalStorage`. It presents a FileLike
    interface to users, but internally (via the :class:`.SharedPath` objects
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
    scratch : os.PathLike
        the scratch directory shared by all objects on this host
    external : :class:`.ExternalStorage`
        external storage resource where objects should eventualy go
    prefix : str
        label for this specific unit
    holding : os.PathLike
        name of the subdirectory of scratch where shared results are
        temporarily stored; default is '.holding'. This must be the same for
        all units within a DAG.
    delete_holding : bool
        whether to delete the contents of the $SCRATCH/$HOLDING/$PREFIX
        directory when this object is deleted
    read_only : bool
        write to prevent NEW files from being written within this shared
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
        read_only: bool = False,
    ):
        self.external = external
        self.scratch = Path(scratch)
        self.prefix = Path(prefix)
        self.read_only = read_only
        self.delete_holding = delete_holding
        self.holding = holding

        self.registry : set[SharedPath] = set()
        # NOTE: the fact that we use $SCRATCH/$HOLDING/$PREFIX instead of
        # $SCRATCH/$PREFIX/$HOLDING is important for 2 reasons:
        # 1. This doesn't take any of the user's namespace from their
        #    $SCRATCH/$PREFIX directory.
        # 2. This allows us to easily use an external FileStorage where the
        #    external storage is exactly the same as this local storage,
        #    meaning that copies to/from the external storage are no-ops.
        #    Use FileStorage(scratch / holding) for that.
        self.shared_dir = self.scratch / holding / prefix
        self.shared_dir.mkdir(exist_ok=True, parents=True)

    def get_other_shared_dir(self, prefix, delete_holding=None):
        """Get a related unit's shared directory.
        """
        if delete_holding is None:
            delete_holding = self.delete_holding

        return SharedRoot(
            scratch=self.scratch,
            external=self.external,
            prefix=prefix,
            holding=self.holding,
            delete_holding=delete_holding,
            read_only=True,
        )

    def transfer_single_file_to_external(self, held_file):
        """Transfer a given file from holding into external storage
        """
        if self.read_only:
            logging.debug("Read-only: Not transfering to external storage")
            return  # early exit

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
        if self.read_only:
            logging.debug("Read-only: Not transfering to external storage")
            return  # early exit

        for obj in self.registry:
            self.transfer_single_file_to_external(obj)

    def cleanup(self):
        """Perform end-of-lifecycle cleanup.
        """
        if self.delete_holding:
            for file in self.registry:
                os.delete(file)
            _delete_empty_dirs(self.shared_dir)

    def register_path(self, shared_path):
        """Register a :class:`.SharedPath` with this :class:`.SharedRoot`.

        This marks a given path as something for this object to manage, by
        loading it into the ``registry``. This way it is tracked such that
        its contents can be transfered to the :class:`.ExternalStorage` and
        such that the local copy can be deleted when it is no longer needed.

        If this objects's :class:`.ExternalStorage` already has data for the
        label associated with the provided :class:`.Sharedpath`, then the
        contents of that should copied to the local path so that it can be
        read by the user.

        Parameters
        ----------
        shared_path: :class:`.SharedPath`
            the path to track
        """
        label_exists = self.external.exists(shared_path.label)

        if self.read_only and not label_exists:
            raise IOError(f"Unable to create '{shared_path.label}'. This "
                          "shared path is read-only.")

        self.registry.add(shared_path)

        # if this is a file that exists, bring it into our subdir
        # NB: this happens even if you're intending to overwrite the path,
        # which is kind of wasteful
        if label_exists:
            scratch_path = self.shared_dir / shared_path.path
            # TODO: switch this to using `get_filename` and `store_path`
            with self.external.load_stream(shared_path.label) as f:
                external_bytes = f.read()
            if scratch_path.exists():
                ... # TODO: something to check that the bytes are the same?
            scratch_path.parent.mkdir(exist_ok=True, parents=True)
            with open(scratch_path, mode='wb') as f:
                f.write(external_bytes)

    def __truediv__(self, path: PathLike):
        return SharedPath(root=self, path=path)

    def __fspath__(self):
        return str(self.shared_dir)

    def __repr__(self):
        return f"SharedRoot({self.scratch}, {self.external}, {self.prefix})"


class SharedPath:
    """PathLike object linking local path with label for external storage.

    On creation, this registers with a :class:`.SharedRoot` that will manage
    the local path and transferring data with its :class:`.ExternalStorage`.
    """
    def __init__(self, root: SharedRoot, path: PathLike):
        self.root = root
        self.path = Path(path)
        self.root.register_path(self)

    def __truediv__(self, path):
        return SharedPath(self.root, self.path / path)

    def __fspath__(self):
        return str(self.root.shared_dir / self.path)

    @property
    def label(self):
        """Label used in :class:`.ExternalStorage` for this path"""
        return str(self.root.prefix / self.path)

    def __repr__(self):
        return f"SharedPath({self.__fspath__()})"

    # TODO: how much of the pathlib.Path interface do we want to wrap?
    # although edge cases may be a pain, we can get most of it with, e.g.:
    # def exists(self): return Path(self).exists()
    # but also, can do pathlib.Path(shared_path) and get hte whole thing
