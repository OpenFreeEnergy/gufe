from contextlib import contextmanager
from pathlib import Path
from typing import Literal

from .externalresource import ExternalStorage


class StorageManager:
    """Manage storage operations for files in a DAG.

    This class provides a context manager for working with storage systems,
    allowing registration, loading, and transfer of files between scratch
    directory and external storage.

    Parameters
    ----------
    scratch_dir : Path
        Path to the scratch directory where files are temporarily stored.
    storage : ExternalStorage
        External storage system for persistent file storage.
    dag_label : str
        Label for the directed acyclic graph (DAG) this storage manager belongs to.
    unit_label : str
        Label for the specific unit within the DAG.

    Attributes
    ----------
    scratch_dir : Path
        Path to the scratch directory.
    storage : ExternalStorage
        External storage system.
    registry : set[str]
        Set of registered filenames to be transferred.
    namespace : str
        Namespace combining dag_label and unit_label for file organization.
    """

    def __init__(
        self,
        scratch_dir: Path,
        storage: ExternalStorage,
        dag_label: str,
        unit_label: str,
    ):
        self.scratch_dir = scratch_dir
        self.storage = storage
        self.registry: set[str] = set()
        self.namespace = f"{dag_label}/{unit_label}"

    @staticmethod
    def append_to_namespace(namespace: str, filename: str) -> str:
        """Append a filenmae to a namespace, mainly used to
        make testing easier.

        Parameters
        ----------
        namespace : str
            The namespace prefix for the file.
        filename : str
            The filename to be appended to the namespace.

        Returns
        -------
        str
            Combined namespace and filename as a storage path.
        """
        # We opt _not_ to use Paths because these aren't actually path objects
        return f"{namespace}/{filename}"

    def register(self, filename: str):
        """Register a filename for later transfer to external storage.

        Parameters
        ----------
        filename : str
            The filename to register for transfer.
        """
        self.registry.add(filename)

    def load(self, filename: str) -> bytes:
        """Load an item from external storage.

        Parameters
        ----------
        filename : str
            The filename to load from external storage.

        Returns
        -------
        bytes
            The content of the loaded file.
        """
        with self.storage.load_stream(filename) as f:
            stored = f.read()
            return stored

    def __contains__(self, filename: str) -> bool:
        return filename in self.registry

    def _transfer(self):
        """Transfer all the files from the files in the internal registry to its
        corresponding :class:`gufe.externalresource.ExternalStorage`.
        """
        for filename in self.registry:
            path = self.scratch_dir / filename
            with open(path, "rb") as f:
                data = f.read()
                self.storage.store_bytes(self.append_to_namespace(self.namespace, filename), data)
