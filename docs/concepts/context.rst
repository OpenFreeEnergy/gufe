Context
=======

``Context`` instances carry the execution environment for individual
``ProtocolUnit`` executions.  They are created by the execution engine just
before a unit's ``_execute`` method is called and then discarded once the unit
returns.  The class lives in ``gufe.protocols.protocolunit`` and acts as a thin
wrapper around two :class:`~gufe.storage.storagemanager.StorageManager`
objects and a scratch directory.


Why Context exists
------------------

``ProtocolUnit`` code frequently needs a few shared facilities:

``scratch``
    A temporary directory that the unit can freely write to while it runs.
    Files written here are considered ephemeral; the engine may delete them as
    soon as the unit finishes.

``shared``
    A :class:`~gufe.storage.storagemanager.StorageManager` backed by
    short–lived :class:`~gufe.storage.externalresource.ExternalStorage`.  Use
    this to hand large files to downstream units without serializing the
    payloads through Python return values.

``permanent``
    Another ``StorageManager`` targeting long–term storage.  Results saved here
    survive beyond the life of the ``ProtocolDAG`` run (for example for
    inspection or for reuse in future extensions).

``stdout`` / ``stderr``
    Optional directories where the engine captures subprocess output triggered
    by the unit.  The directories are removed automatically when the context
    closes.

Keeping these handles bundled together and managed by a context manager lets
``ProtocolUnit`` implementers focus on domain logic while the engine ensures
storage gets flushed and temporary directories are cleaned.


Lifecycle
---------

``Context`` implements ``__enter__``/``__exit__`` so the execution engine can
rely on ``with Context(...) as ctx`` semantics internally.  When ``__exit__``
triggers the ``shared`` and ``permanent`` storage managers flush pending
transfers back to their underlying ``ExternalStorage`` instances.  Any
``stdout`` or ``stderr`` capture directories are also removed.

This means a ``ProtocolUnit`` should treat ``ctx.shared`` and
``ctx.permanent`` as handles that stay valid only for the duration of
``_execute``; once the method returns, the engine is responsible for moving the
files into their durable homes.


Using Context inside ProtocolUnits
----------------------------------

Every ``ProtocolUnit._execute`` definition must accept ``ctx: Context`` as its
first argument.  Typical usage looks like the example below.

.. code-block:: python

    from gufe import ProtocolUnit, Context

    class SimulationUnit(ProtocolUnit):

        @staticmethod
        def _execute(ctx: Context, *, setup_result, lambda_window, settings):
            scratch_path = ctx.scratch / f"lambda_{lambda_window}"
            scratch_path.mkdir(exist_ok=True)

            # Read upstream artifacts from ctx.shared
            system_file = setup_result.outputs["system_file"]
            topology_file = setup_result.outputs["topology_file"]

            # Produce large payloads into ctx.shared so downstream units can
            # access them without serialization
            result_path = ctx.shared / f"window_{lambda_window}.npz"
            save_simulation_outputs(result_path)

            # Return only lightweight metadata
            return {
                "lambda_window": lambda_window,
                "result_path": str(result_path),
            }


Choosing between shared and permanent storage
---------------------------------------------

Both ``ctx.shared`` and ``ctx.permanent`` expose the same ``StorageManager``
API but they serve different audiences:

``ctx.shared``
    Optimized for communication between units in the same DAG execution.  The
    execution backend is free to prune these assets once no downstream unit
    references them.

``ctx.permanent``
    Intended for outputs that should survive beyond the immediate DAG, such as
    user-facing reports or artifacts that will seed future ``extends`` runs.

As a rule of thumb, prefer ``ctx.shared`` unless you have a clear requirement
to keep the data after the ``Protocol`` run concludes.  Small scalar values or
lightweight metadata should still be returned directly from ``_execute`` so
they become part of the ``ProtocolUnitResult`` record.


Interaction with Protocols
--------------------------

``Protocol`` instances do not instantiate ``Context`` directly; they declare
``ProtocolUnit`` objects via ``Protocol._create``.  When an execution backend
(for example, ``gufe-client`` or a workflow manager) walks the resulting
``ProtocolDAG`` it constructs a ``Context`` for each unit using the DAG label,
unit label, scratch directory, and the configured
``ExternalStorage`` implementations.  The backend might provide different
``ExternalStorage`` implementations (e.g., local filesystem, object store,
cluster scratch) depending on where the work runs, but the ``Context`` API seen
by ``ProtocolUnit`` authors stays consistent.

Because the execution backend is in charge of creating the contexts, protocol
authors can rely on ``ctx`` always being populated with valid storage managers
and paths that are safe to write to from distributed workers.


Practical tips
--------------

* Keep return values small.  Record numeric summaries, filenames, and status in
  the ``dict`` returned by ``_execute``; stream heavy data through
  ``ctx.shared`` or ``ctx.permanent``.
* Clean up unit-specific scratch files as you go if they consume lots of disk;
  the execution engine only guarantees cleanup once ``_execute`` finishes.
* If you launch subprocesses, direct their output into ``ctx.stdout`` and
  ``ctx.stderr`` (when provided) to make debugging easier.
* ``StorageManager`` exposes convenience helpers such as
  :meth:`~gufe.storage.storagemanager.StorageManager.store_path` for shuttling
  directories or streams; use these instead of manual ``shutil`` calls.

By following these patterns, ``ProtocolUnit`` implementations remain portable
across execution backends while benefiting from consistent lifecycle and
storage management.


Migrating from legacy Context usage
-----------------------------------

Before ``Context`` was rewritten in :mod:`gufe 0.15`, it was a simple data
class with two ``pathlib.Path`` handles: ``scratch`` and ``shared``.  Existing
protocols adopted a variety of implicit conventions around those attributes.
The current implementation keeps backwards compatibility at the interface
level—``ctx.scratch`` and ``ctx.shared`` still exist—but their behavior is more
structured because they are ``StorageManager`` instances backed by
``ExternalStorage``.

Follow this checklist when migrating old protocols:

1. **Swap file paths to StorageManager APIs.**  Calls like ``ctx.shared /
   "filename"`` should be replaced with the helper methods offered by
   :class:`StorageManager`.  For example::

       path = ctx.shared.scratch_dir / "myfile.dat"
       ctx.shared.register("myfile.dat")

   or use :meth:`StorageManager.store_path` when moving directories.  Register
   every artifact that must be transferred off the worker; unregistered files
   stay local and will be deleted.

2. **Avoid storing heavy objects in Python outputs.**  Older protocols often
   returned raw ``Path`` objects pointing at scratch files.  Instead, register
   the file with the storage manager and return the storage key (a string) from
   ``_execute`` as shown in the example above.  Downstream units can then call
   ``ctx.shared.storage.load_stream`` or ``StorageManager.load``.

3. **Handle ``ctx.permanent``.**  There was no equivalent in the legacy API.
   Decide which results must persist between DAG executions and write them via
   ``ctx.permanent``.  For migration you can start by mirroring whatever used
   to live in ``ctx.shared`` and refine later.

4. **Expect automatic cleanup.**  Old contexts typically left stdout/stderr
   directories around.  The new context removes these when the unit finishes.
   If your code tried to re-read the capture directories after ``_execute``
   returned, move that logic earlier or rely on the logged data captured in the
   ``ProtocolUnitResult``.

5. **Stop constructing Context manually.**  Some bespoke execution scripts once
   instantiated ``Context(scratch=..., shared=...)`` by hand.  That pattern is
   obsolete because the constructor now requires ``ExternalStorage`` objects.
   Instead, rely on the execution backend (``gufe-client``, custom schedulers,
   etc.) to build contexts.  For unit tests use the helpers in
   ``gufe.storage.externalresource`` (e.g., ``MemoryStorage``) to create the
   necessary storage instances.

During the transition it is safe to read ``ctx.shared`` as a ``Path`` so long as
you treat it as read-only and only for temporary files; the attribute still
implements ``pathlib.Path`` thanks to the embedded ``StorageManager``'s
``scratch_dir``.  Prefer using the explicit ``scratch_dir`` attribute to make
the intent clear when touching raw filesystem paths.
