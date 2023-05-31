from os import PathLike
from pathlib import Path
from contextlib import contextmanager

from .externalstorage import ExternalStorage, FileStorage

def _storage_path_conflict(external, path):
    # this is a little brittle; I don't like hard-coding the class here
    if isinstance(external, FileStorage):
        root = Path(external.root_dir)
    else:
        return False

    return False


class _DAGStorageManager:
    """Context manager to handle details of storage lifecycle.

    Making this a separate class ensures that ``running_unit`` is always
    called within the context of a given DAG. This is usually not created
    directly; instead, it is created (and used) with its ``running_dag``
    classmethod, typically from within a ``StorageManager``.
    """
    def __init__(self, storage_manager, dag_label):
        self.manager = storage_manager
        self.dag_label = dag_label
        self.permanents = []

    @classmethod  # NB: classmethod must be on top
    @contextmanager
    def running_dag(cls, storage_manager, dag_label):
        """DAG level of the storage lifecycle

        When the DAG is completed, transfer everything to the permanent
        storage, and delete the holding area for permanent (if we are
        supposed to).

        This is not usually called by users; instead it is called from
        within the ``StorageManager``.
        """
        dag_manager = cls(storage_manager, dag_label)
        try:
            yield dag_manager
        finally:
            for permanent in dag_manager.permanents:
                permanent.transfer_holding_to_external()

            if not dag_manager.manager.keep_holding:
                for d in dag_manager.permanents:
                    d.cleanup()

    @contextmanager
    def running_unit(self, unit):
        """Unit level of the storage lifecycle.

        This provides the holding directories used for scratch, shared, and
        permanent. At the end of the unit, it transfers anything from shared
        to the real shared external storage, cleans up the scratch
        directory and the shared holding directory.
        """
        unit_label = unit.key
        scratch = self.manager.get_scratch(self.dag_label, unit_label)
        shared = self.manager.get_shared(self.dag_label, unit_label)
        permanent = self.manager.get_permanent(self.dag_label, unit_label)
        try:
            yield scratch, shared, permanent
        finally:
            # TODO: should some of this be in an else clause instead?
            self.permanents.append(permanent)
            shared.transfer_holding_to_external()
            # everything in permanent must also be available in shared
            for file in permanent.registry:
                shared.transfer_single_file_to_external(file)
            scratch_conflict = _storage_path_conflict(shared.external,
                                                      scratch)
            if not self.manager.keep_scratch and not scratch_conflict:
                shutil.rmtree(scratch)

            shared_conflict = _storage_path_conflict(shared.external,
                                                     shared)

            if not self.manager.keep_holding:
                shared.cleanup()


class StorageManager:
    """Tool to manage the storage lifecycle during a DAG.

    This object primarily contains the logic for getting the holding
    directories. A separate class, in the ``DAGContextClass`` variable,
    handles the logic for the context managers.
    """
    def __init__(
        self,
        scratch_root: PathLike,
        shared_root: ExternalStorage,
        permanent_root: ExternalStorage,
        *,
        keep_scratch: bool = False,
        keep_holding: bool = False,
        holding: PathLike = Path(".holding"),
        DAGContextClass: type = _DAGStorageManager,
    ):
        self.scratch_root = scratch_root
        self.shared_root = shared_root
        self.permanent_root = permanent_root
        self.keep_scratch = keep_scratch
        self.keep_holding = keep_holding
        self.holding = holding
        self.DAGContextClass = DAGContextClass

    def get_scratch(self, dag_label: str , unit_label: str) -> Path:
        """Get the path for this unit's scratch directory"""
        scratch = self.scratch_root / dag_label / "scratch" / unit_label
        scratch.mkdir(parents=True, exist_ok=True)
        return scratch

    def get_permanent(self, dag_label, unit_label):
        """Get the object for this unit's permanent holding directory"""
        return SharedRoot(
            scratch=self.scratch_root / dag_label,
            external=self.permanent_root,
            prefix=unit_label,
        )

    def get_shared(self, dag_label, unit_label):
        """Get the object for this unit's shared holding directory"""
        return SharedRoot(
            scratch=self.scratch_root / dag_label,
            external=self.shared_root,
            prefix=unit_label
        )

    def running_dag(self, dag_label: str):
        """Return a context manager that handles storage.

        For simple use cases, this is the only method a user needs to call.
        Usage is something like:

        .. code::

            with manager.running_dag(dag_label) as dag_ctx:
                for unit in dag_ordered_units:
                    with dag_ctx.running_unit(unit) as dirs:
                        scratch, shared, permanent = dirs
                        # run the unit
        """
        return self.DAGContextClass.running_dag(self, dag_label)
