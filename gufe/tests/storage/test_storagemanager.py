# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import pathlib
import tempfile
from unittest.mock import MagicMock, patch

import pytest
from pluggy import _manager

from gufe.storage.externalresource import FileStorage, MemoryStorage
from gufe.storage.externalresource.base import ExternalStorage
from gufe.storage.storagemanager import StorageManager


@pytest.fixture
def tmp_scratch_dir():
    """Create a temporary scratch directory."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield pathlib.Path(tmp_dir)


@pytest.fixture
def memory_storage():
    """Create a MemoryStorage instance for testing."""
    return MemoryStorage()


@pytest.fixture
def file_storage(tmp_path):
    """Create a FileStorage instance for testing."""
    return FileStorage(tmp_path)


@pytest.fixture
def storage_manager(tmp_scratch_dir, memory_storage):
    """Create a StorageManager instance with MemoryStorage."""
    return StorageManager(tmp_scratch_dir, memory_storage, dag_label="MEM", unit_label="1")


@pytest.fixture
def storage_manager_file_storage(tmp_scratch_dir, file_storage):
    """Create a StorageManager instance with FileStorage."""
    return StorageManager(tmp_scratch_dir, file_storage, dag_label="FILE", unit_label="1")


class TestStorageManager:
    """Test the StorageManager class."""

    def test_init(self, tmp_scratch_dir, memory_storage):
        """Test StorageManager initialization."""
        manager = StorageManager(tmp_scratch_dir, memory_storage, dag_label="MEM", unit_label="1")

        assert manager.scratch_dir == tmp_scratch_dir
        assert manager.storage == memory_storage
        assert isinstance(manager.registry, set)
        assert len(manager.registry) == 0
        assert manager.namespace == "MEM/1"

    def test_init_with_file_storage(self, tmp_scratch_dir, file_storage):
        """Test StorageManager initialization with FileStorage."""
        manager = StorageManager(tmp_scratch_dir, file_storage, dag_label="FILE", unit_label="1")

        assert manager.scratch_dir == tmp_scratch_dir
        assert manager.storage == file_storage
        assert isinstance(manager.registry, set)

    def test_convert_to_namespace(self):
        filename = "test_file.txt"

        out = StorageManager.append_to_namespace("MEM/1", filename)
        assert out == "MEM/1/test_file.txt"

    def test_register(self, storage_manager):
        """Test registering a filename."""
        filename = "test_file.txt"

        # Initially registry should be empty
        assert filename not in storage_manager.registry
        assert filename not in storage_manager

        # Register the file
        out = storage_manager.register(filename)
        assert out == storage_manager.append_to_namespace(storage_manager.namespace, filename)

        # Check it's now in the registry
        assert filename in storage_manager.registry
        assert filename in storage_manager

    def test_register_multiple_files(self, storage_manager):
        """Test registering multiple filenames."""
        files = ["file1.txt", "file2.txt", "dir/file3.txt"]

        for filename in files:
            storage_manager.register(filename)

        # Check all files are registered
        for filename in files:
            assert filename in storage_manager.registry
            assert filename in storage_manager

        # Check registry size
        assert len(storage_manager.registry) == 3

    def test_register_duplicate_file(self, storage_manager):
        """Test registering the same file multiple times."""
        filename = "duplicate.txt"

        # Register twice
        storage_manager.register(filename)
        storage_manager.register(filename)

        # Should only appear once (set behavior)
        assert filename in storage_manager.registry
        assert len(storage_manager.registry) == 1

    def test_contains(self, storage_manager):
        """Test the __contains__ method."""
        filename = "contains_test.txt"

        # Initially should not contain
        assert filename not in storage_manager

        # Register and check again
        storage_manager.register(filename)
        assert filename in storage_manager

        # Check a non-existent file
        assert "nonexistent.txt" not in storage_manager

    def test_transfer_empty_registry(self, storage_manager):
        """Test _transfer with empty registry."""
        # Should not raise any errors
        storage_manager._transfer()

        # Storage should remain empty
        assert list(storage_manager.storage) == []

    def test_transfer_with_files(self, storage_manager, tmp_scratch_dir):
        """Test _transfer with actual files."""
        # Create test files in scratch directory
        test_files = {"test1.txt": b"Hello World", "test2.txt": b"Another file", "subdir/test3.txt": b"Nested file"}

        # Create files and register them
        for filename, content in test_files.items():
            file_path = tmp_scratch_dir / filename
            file_path.parent.mkdir(parents=True, exist_ok=True)
            file_path.write_bytes(content)
            storage_manager.register(filename)

        # Transfer files
        storage_manager._transfer()

        # Check files are now in storage
        for filename in test_files:
            namespaced_filename = StorageManager.append_to_namespace(storage_manager.namespace, filename)
            assert storage_manager.storage.exists(namespaced_filename)

            # Verify content
            with storage_manager.storage.load_stream(namespaced_filename) as f:
                stored_content = f.read()
            assert stored_content == test_files[filename]

    def test_transfer_with_file_storage(self, storage_manager_file_storage, tmp_scratch_dir):
        """Test _transfer with FileStorage."""
        # Create test file
        filename = "transfer_test.txt"
        content = b"Test content for file storage"

        # Create file in scratch directory
        file_path = tmp_scratch_dir / filename
        file_path.write_bytes(content)
        storage_manager_file_storage.register(filename)

        # Transfer file
        storage_manager_file_storage._transfer()

        # Check file exists in file storage
        namespaced_filename = StorageManager.append_to_namespace(storage_manager_file_storage.namespace, filename)
        assert storage_manager_file_storage.storage.exists(namespaced_filename)

        # Verify content
        expected_path = storage_manager_file_storage.storage.root_dir / namespaced_filename
        assert expected_path.exists()
        assert expected_path.read_bytes() == content

    def test_transfer_file_not_found(self, storage_manager):
        """Test _transfer when registered file doesn't exist."""
        filename = "nonexistent.txt"

        # Register non-existent file
        storage_manager.register(filename)

        # Should raise FileNotFoundError when trying to transfer
        with pytest.raises(FileNotFoundError):
            storage_manager._transfer()

    def test_registry_persistence(self, storage_manager, tmp_scratch_dir):
        """Test that registry persists across operations."""
        filename = "persistent.txt"

        # Create the file first
        file_path = tmp_scratch_dir / filename
        file_path.write_bytes(b"test content")

        # Register file
        storage_manager.register(filename)
        assert filename in storage_manager.registry

        # Perform other operations
        storage_manager._transfer()  # Should not clear registry

        # Registry should still contain the file
        assert filename in storage_manager.registry

    def test_load(self, storage_manager, tmp_scratch_dir):
        filename = "item.txt"
        file_path = tmp_scratch_dir / filename
        content = b"test content"
        file_path.write_bytes(content)

        # Register file
        storage_manager.register(filename)
        assert filename in storage_manager.registry

        # Perform other operations
        storage_manager._transfer()  # Should not clear registry

        # Registry should still contain the file
        assert filename in storage_manager.registry
        namespaced_filename = StorageManager.append_to_namespace(storage_manager.namespace, filename)
        out = storage_manager.load(namespaced_filename)
        assert out == content

    def test_prepopulated_storage_load(self, file_storage: ExternalStorage, tmp_scratch_dir):
        # This test highlights a case where a unit wants to load something
        # from shared storage. We basically load content into the medium that we want to fetch later.

        # Store content
        contents = b"test content"
        name = "test.txt"
        file_storage.store_bytes(name, contents)
        # Provide this prepopulated storage to a manager
        manager = StorageManager(tmp_scratch_dir, storage=file_storage, dag_label="TEST", unit_label="1")
        # Validate that content is loaded
        out = manager.load(name)
        assert out == contents

        # Validate that the registry has not loaded the item
        assert len(manager.registry) == 0
