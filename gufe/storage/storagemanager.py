from __future__ import annotations
from os import PathLike
from pathlib import Path
from contextlib import contextmanager
import shutil

from gufe.utils import delete_empty_dirs

from typing import Type

from .externalresource import ExternalStorage, FileStorage
from .stagingregistry import SharedStaging, PermanentStaging
from .stagingregistry import StagingPath  # typing


class StorageManager:
    def __init__(
        self,
        scratch_root: PathLike,
        shared_root: ExternalStorage,
        permanent_root: ExternalStorage,
        *,
        keep_scratch: bool = False,
        keep_staging: bool = False,
        keep_shared: bool = False,
        keep_empty_dirs: bool = False,
        staging: PathLike = Path(".staging"),
    ):
        self.scratch_root = Path(scratch_root)
        self.shared_root = shared_root
        self.permanent_root = permanent_root
        self.keep_scratch = keep_scratch
        self.keep_staging = keep_staging
        self.keep_shared = keep_shared
        self.staging = staging
        self.keep_empty_dirs = keep_empty_dirs

        # these are used to track what files can be deleted from shared if
        # keep_shared is False
        self.shared_xfer: set[StagingPath] = set()
        self.permanent_xfer: set[StagingPath] = set()

        self.permanent_staging = PermanentStaging(
            scratch=self.scratch_root,
            external=self.permanent_root,
            shared=self.shared_root,
            staging=self.staging,
            keep_empty_dirs=keep_empty_dirs,
        )

        self.shared_staging = SharedStaging(
            scratch=self.scratch_root,
            external=self.shared_root,
            staging=self.staging,
            keep_empty_dirs=keep_empty_dirs,
        )

    def make_label(self, dag_label, unit_label, attempt, **kwargs):
        """

        The specific executor may change this by making a very simple
        adapter subclass and overriding this method, which can take
        arbitrary additional kwargs that may tie it to a specific executor.
        """
        return f"{dag_label}/{unit_label}_attempt_{attempt}"

    @property
    def _scratch_base(self):
        return self.scratch_root / "scratch"

    def _scratch_loc(self, dag_label, unit_label, attempt, **kwargs):
        label = self.make_label(dag_label, unit_label, attempt)
        return self._scratch_base / label

    @contextmanager
    def running_dag(self, dag_label):
        # TODO: remove (or use) dag_label
        try:
            yield self
        finally:
            # import pdb; pdb.set_trace()
            # clean up after DAG completes
            self.permanent_staging.transfer_staging_to_external()

            if not self.keep_staging:
                self.permanent_staging.cleanup()

            if not self.keep_shared:
                # we'd like to do something like loop over
                # self.shared_xfer - self.permanent_xfer; however,
                # StagedPaths have different staging registries. This gives
                # the set of paths we do want to delete
                perm_xfer_paths = {p.as_path() for p in self.permanent_xfer}
                shared_xfer_to_delete = {
                    p for p in self.shared_xfer
                    if p.as_path() not in perm_xfer_paths
                }

                for file in shared_xfer_to_delete:
                    self.shared_root.delete(file.label)

                for file in self.permanent_xfer:
                    if self.shared_root != self.permanent_root:
                        self.shared_root.delete(file.label)

            if not self.keep_empty_dirs:
                delete_empty_dirs(self._scratch_base, delete_root=False)

    @contextmanager
    def running_unit(self, dag_label, unit_label, **kwargs):
        scratch = self._scratch_loc(dag_label, unit_label, **kwargs)
        label = self.make_label(dag_label, unit_label, **kwargs)
        scratch.mkdir(parents=True, exist_ok=True)
        shared = self.shared_staging / label
        permanent = self.permanent_staging / label
        try:
            yield scratch, shared, permanent
        finally:
            # import pdb; pdb.set_trace()
            # clean up after unit

            # track the files that were in shared so that we can delete them
            # at the end of the DAG if requires
            shared_xfers = self.shared_staging.transfer_staging_to_external()
            self.shared_xfer.update(set(shared_xfers))

            # everything in permanent should also be in shared
            for file in self.permanent_staging.registry:
                self.shared_staging.transfer_single_file_to_external(file)
                self.permanent_xfer.add(file)

            if not self.keep_scratch:
                shutil.rmtree(scratch)

            if not self.keep_staging:
                self.shared_staging.cleanup()
