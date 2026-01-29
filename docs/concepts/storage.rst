.. _concepts-storage:

How storage is handled in **gufe**
==================================

**gufe** abstracts storage into a reusable interface using the :class:`.ExternalStorage` abstract base class.
This abstraction enables the storage of any file or byte stream using various storage backends without changing application code.

Overview
--------

The storage system is designed to handle files (or byte data) that need to be stored in some location.
Instead of embedding the data, objects can store a reference (a unique string indicating the object's location, such as a path) to where the data is stored externally. This approach provides several benefits:

* **Efficiency**: Large objects don't need to be serialized multiple times
* **Flexibility**: Different storage backends (local filesystem, cloud storage, in-memory) can be used interchangeably
* **Deduplication**: The same data can be referenced by multiple objects
* **Lazy Loading**: Data is only loaded when needed

The Storage Architecture
-------------------------

The storage system consists of several key components:

``ExternalStorage`` Base Class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :class:`.ExternalStorage` abstract base class defines the interface that all storage implementations must provide. This class provides:

* **Store operations**: ``store_bytes()`` and ``store_path()`` to store data
* **Load operations**: ``load_stream()`` to retrieve data as a stream
* **Management**: ``exists()``, ``delete()``, and ``iter_contents()`` for managing stored data

All storage operations use a *location* string as an identifier for the stored data.

Storage Implementations
-----------------------

**gufe** provides several built-in implementations of :class:`.ExternalStorage`:

``FileStorage``
~~~~~~~~~~~~~~~

The :class:`.FileStorage` implementation stores data on the local filesystem. It requires a root directory path and organizes stored files using the location string as a relative path:

.. code-block:: python

    from pathlib import Path
    from gufe.storage.externalresource import FileStorage

    # Create a file storage backend
    storage = FileStorage(root_dir=Path("/path/to/storage"))

    # Store some data
    data = b"Hello, World!"
    storage.store_bytes("datasets/sample1.txt", data)

    # Check if data exists
    if storage.exists("datasets/sample1.txt"):
        # Load the data
        with storage.load_stream("datasets/sample1.txt") as stream:
            loaded_data = stream.read()
            assert loaded_data == data

    # Delete the data
    storage.delete("datasets/sample1.txt")

``FileStorage`` automatically creates any necessary parent directories when storing files.

``MemoryStorage``
~~~~~~~~~~~~~~~~~

The :class:`.MemoryStorage` implementation stores data in a Python dictionary. This is primarily useful for testing and prototyping:

.. code-block:: python

    from gufe.storage.externalresource import MemoryStorage

    # Create an in-memory storage backend
    storage = MemoryStorage()

    # Store some data
    data = b"Hello, World!"
    storage.store_bytes("datasets/sample1.txt", data)

    # Load the data back
    with storage.load_stream("datasets/sample1.txt") as stream:
        loaded_data = stream.read()

.. warning::
    ``MemoryStorage`` is not intended for production use and all data is lost when the Python process exits.

Implementing Custom Storage Backends
-------------------------------------

To create a custom storage backend, subclass :class:`.ExternalStorage` and implement all the abstract methods:

.. code-block:: python

    from gufe.storage.externalresource.base import ExternalStorage
    from typing import ContextManager

    class MyCustomStorage(ExternalStorage):
        """A custom storage implementation."""

        def _store_bytes(self, location: str, byte_data: bytes):
            """Store bytes at the given location."""
            # Implement storage logic
            pass

        def _store_path(self, location: str, path):
            """Store a file at the given path."""
            # Implement storage logic
            pass

        def _load_stream(self, location: str) -> ContextManager:
            """Return a context manager that yields a bytes-like object."""
            # Implement loading logic
            pass

        def _exists(self, location: str) -> bool:
            """Check if data exists at the location."""
            # Implement existence check
            pass

        def _delete(self, location: str):
            """Delete data at the location."""
            # Implement deletion logic
            pass

        def _get_filename(self, location: str) -> str:
            """Return a filename for the location."""
            # Implement filename generation
            pass

        def _iter_contents(self, prefix: str = ""):
            """Iterate over stored locations matching the prefix."""
            # Implement iteration logic
            pass

        def _get_hexdigest(self, location: str) -> str:
            """Return MD5 hexdigest of the data (optional override)."""
            # Can override for performance improvements
            pass

.. note::
    All storage methods should be blocking operations, even if the underlying storage backend supports asynchronous operations.

StorageManager
--------------

The :class:`.StorageManager` class provides a higher-level interface for managing storage operations within a computational workflow.
.. note::

    ``StorageManager`` is largely used by the :class:`.Context` class and should not be instantiated in protocols.
    In general, protocol developers will only use the ``register`` and ``load`` functions.

It handles the transfer of files between a scratch directory and external storage (such as shared or permenant storage):

.. code-block:: python

    from pathlib import Path
    from gufe.storage import StorageManager
    from gufe.storage.externalresource import FileStorage

    # Set up storage
    storage = FileStorage(root_dir=Path("/path/to/storage"))
    scratch_dir = Path("/path/to/scratch")

    # Create a storage manager for a specific DAG and unit
    manager = StorageManager(
        scratch_dir=scratch_dir,
        storage=storage,
        dag_label="my_experiment",
        unit_label="transformation_1"
    )

    # Register files for later transfer
    out = manager.register("trajectory.dcd")
    out2 = manager.register("results.json")
    # Note: out and out2 are pre-namespaced values that allow storage items to be passed around

    # Transfer all registered files to external storage
    manager._transfer()

    # Load files from external storage
    trajectory_data = manager.load(out)
    results_json = manager.load(out2)

The ``StorageManager`` uses a namespace combining the ``dag_label`` and ``unit_label`` to organize files in the external storage backend.
To see how these work in practice see our documentation on :doc:`context`.


Error Handling
--------------

The storage system defines several exceptions in :mod:`gufe.storage.errors`:

* :class:`.ExternalResourceError`: Base class for storage-related errors
* :class:`.MissingExternalResourceError`: Raised when attempting to access non-existent data
* :class:`.ChangedExternalResourceError`: Raised when metadata verification fails

These exceptions can be caught and handled appropriately in application code:

.. code-block:: python

    from gufe.storage.errors import MissingExternalResourceError

    try:
        with storage.load_stream("nonexistent_file.txt") as stream:
            data = stream.read()
    except MissingExternalResourceError:
        print("File not found in storage")
