from contextlib import contextmanager
from pathlib import Path
from typing import Literal

from .externalresource import ExternalStorage


class StorageManager:
    """This class exists to be context manager for working with storage systems."""

    def __init__(self, scratch_path: Path, storage: ExternalStorage):
        self.scratch_path = scratch_path
        self.storage = storage
        self.registry: set[str] = set()

    @property
    def encoder(self):
        """The method used for how to encode information to storage."""
        raise NotImplemented

    @property
    def decoder(self):
        """The method used for how to decode information to storage."""
        raise NotImplemented

    def _register(self, filename: str):
        """Register a filename to a given store so it can be moved later."""
        # TODO: Check if the handler already exists.
        self.registry.add(filename)

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
                self.storage.store_bytes(filename, data)

    @contextmanager
    def running_unit(self, dag_label, unit_label):
        # TODO: This needs to yield the correct gufe Context object
        raise NotImplemented
