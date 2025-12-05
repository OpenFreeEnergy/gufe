from contextlib import contextmanager
from pathlib import Path
from typing import Literal

from .externalresource import ExternalStorage


class StorageManager:
    """This class exists to be context manager for working with storage systems."""

    def __init__(
        self,
        scratch_path: Path,
        storage: ExternalStorage,
        dag_label: str,
        unit_label: str,
    ):
        self.scratch_path = scratch_path
        self.storage = storage
        self.registry: set[str] = set()
        self.namespace = f"{dag_label}/{unit_label}"

    def _convert_to_namespace(self, filename: str) -> str:
        # We opt _not_ to use Paths because these aren't actually path objects
        return f"{self.namespace}/{filename}"

    def register(self, filename: str):
        """Register a filename to a given store so it can be moved later."""
        self.registry.add(filename)

    def load(self, filename: str) -> bytes:
        """Load an item from external storage"""
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
            path = self.scratch_path / filename
            with open(path, "rb") as f:
                data = f.read()
                self.storage.store_bytes(self._convert_to_namespace(filename), data)

    @contextmanager
    def running_unit(self, dag_label, unit_label):
        # TODO: This needs to yield the correct gufe Context object
        raise NotImplemented
