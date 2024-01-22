from gufe.storage.storagemanager import StorageManager
from gufe.storage.stagingserialization import StagingPathSerialization
from gufe.storage.externalresource import ExternalStorage
from contextlib import contextmanager
from pathlib import Path
from os import PathLike

from gufe.protocols.protocoldag import Context

class PerUnitStorageManager(StorageManager):
    """Variant to use when doing only a single process per unit"""
    def __init__(
        self,
        scratch_root: PathLike,
        shared_root: ExternalStorage,
        permanent_root: ExternalStorage,
        *,
        keep_scratch: bool = False,
        keep_staging: bool = False,
        keep_shared: bool = False,
        staging: PathLike = Path(".staging"),
        delete_empty_dirs: bool = True,
    ):
        super().__init__(
            scratch_root=scratch_root,
            shared_root=shared_root,
            permanent_root=permanent_root,
            keep_scratch=keep_scratch,
            keep_staging=keep_staging,
            keep_shared=keep_shared,
            delete_empty_dirs=delete_empty_dirs,
        )
        # TODO: move this to the base class
        self.serialization = StagingPathSerialization(self)

    @property
    def json_encoder(self):
        self.serialization.refresh_handler()
        return self.serialization.encoder

    @property
    def json_decoder(self):
        self.serialization.refresh_handler()
        return self.serialization.decoder

    @contextmanager
    def running_dag(self, dag_label):
        yield self

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
            self.shared_staging.transfer_staging_to_external()
            for file in self.permanent_staging.registry:
                self.shared_staging.transfer_single_file_to_external(file)

            self.permanent_staging.transfer_staging_to_external()
