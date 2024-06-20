The GUFE Storage System
=======================

Storage lifetimes
-----------------

The storage system in GUFE is heavily tied to the GUFE protocol system. The
key concept here is that the different levels of the GUFE protocol system;
campaign, DAG, and unit; inherently imply different lifetimes for the data
that is created. Those different lifetimes define the stages of the GUFE
storage system.

In an abstract sense, as used by protocol developers, these three levels
correspond to three lifetimes of data:

* ``scratch``: This is temporary data that is only needed for the lifetime
  of a :class:`.ProtocolUnit`. This data is not guaranteed to be available
  beyond the single :class:`.ProtocolUnit` where it is created, but may be
  reused within that :class:`.ProtocolUnit`.
* ``shared``: This is data that is shared between different units in a
  :class:`.ProtocolDAG`. For example, a single equilibration stage might be
  shared between multiple production runs. The output snapshot of the
  equilibration would be suitable for as something to put in ``shared``
  data. This data is guaranteed to be present from when it is created until
  the end of the :class:`.ProtocolDAG`, but is not guaranteed to exist after
  the :class:`.ProtocolDAG` terminates.
* ``permanent``: This is the data that will be needed beyond the scope of a
  single rough estimate of the calculation. This could include anything that
  an extension of the simulation would require, or things that require
  network-scale analysis. Anything stored here will be usable after the
  calculation has completed.

The ``scratch`` area is always a local directory. However, ``shared`` and
``permanent`` can be external (remote) resources, using the
:class:`.ExternalResource` API.

As a practical matter, the GUFE storage system can be handled with a
:class:`.StorageManager`. This automates some aspects of the transfer
between stages of the GUFE storage system, and simplifies the API for
protocol authors.  In detail, this provides protocol authors with
``PathLike`` objects for ``scratch``, ``shared``, and ``permanent``. All
three of these objects actually point to special subdirectories of the
local scratch space for a specific unit, but are managed by context
managers at the executor level, which handle the process of moving objects
from local staging directories to the actual ``shared`` and ``permanent``
locations, which can be external resources.


External resource utilities
---------------------------

For flexible data storage, GUFE defines the :class:`.ExternalResource` API,
which allows data be stored/loaded in a way that is agnostic to the
underlying data store, as long as the store can be represented as a
key-value store. Withing GUFE, we provide several examples, including
:class:`.FileStorage` and :class:`.MemoryStorage` (primarily useful for
testing.) The specific ``shared`` and ``permanent`` resources, as provided
to the executor, can be instances of an :class:`.ExternalResource`.

.. note::

   The ``shared`` space must be a resource where an uploaded object is
   instantaneously available, otherwise later protocol units may fail if the
   shared result is unavailable. This means that files or approaches based
   on ``scp`` or ``sftp`` are fine, but things like cloud storage, where the
   existence of a new document may take time to propagate through the
   network, are not recommended for ``shared``.


Details: Manangement of the Storage Lifetime
--------------------------------------------

The concepts of the storage lifetimes are important for protocol authors,
but details of implementation are left to the specific executor. In order to
facilitate correct treatment of the storage lifecycle, GUFE provides a few
helpers. The information in this section is mostly of interest to authors of
executors. The helpers are:

* :class:`.StorageManager`: This is the overall façade interface for
  interacting with the rest of the storage lifecycle tools. It provides two
  methods to generate context managers; one for the :class:`.ProtocolDAG`
  level of the lifecycle, and one for the :class:`.ProtocoUnit` level of the
  lifecycle. This class is designed for the use case that the entire DAG is
  run in serial within a single process. Subclasses of this can be created
  for other execution architectures, where the main logic changes would be
  in the methods that return those context managers.
* :class:`.StagingRegistry`: This handles the logic around staging paths
  within a :class:`.ProtocolUnit`. Think of this as an abstract
  representation of a local directory. Paths within it register with it, and
  it handles deletion of the temporary local files when not needed, as well
  as the download of remote files when necessary for reading. There are two
  important subclasses of this: :class:`.SharedStaging` for a ``shared``
  resource, and :class:`.PermanentStaging` for a ``permanent`` resource.
* :class:`.StagingPath`: This represents a file within the
  :class:`.StagingRegistry`. It contains both the key (label) used in the
  key-value store, as well as the actual local path to the file. When its
  ``__fspath__`` method is called, it registers itself with its
  :class:`.StagingRegistry`, which handles managing it over its lifecycle.

In practice, the executor uses the :class:`.StorageManager` to create a
:class:`.DAGContextManager` at the level of a DAG, and then uses the
:class:`.DAGContextManager` to create a context to run a unit. That context
creates a :class:`.SharedStaging` and a :class:`.PermanentStaging`
associated with the specific unit. Those staging directories, with the
scratch directory, are provided to the :class:`.ProtocolUnit`, so that
these are the only objects protocol authors need to interact with.
