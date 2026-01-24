``ProtocolUnit`` execution ``Context``
=======

:class:`.Context` instances carry the execution environment for individual :class:`.ProtocolUnit` executions.
They are created by the execution engine just before a unit is excuted and discarded once the unit returns.
The class acts as a thin wrapper around two :class:`.StorageManager` objects (shared and permanent) and a scratch directory.


Why Context exists
------------------

``ProtocolUnit`` code frequently needs a few shared facilities:

``scratch``
    A temporary directory that the unit can freely write to while it runs.
    Files written here are considered ephemeral; the engine may delete them as
    soon as the unit finishes.

``shared``
    A :class:`.StorageManager` backed by
    short–lived :class:`.ExternalStorage`. Use
    this to hand large files to downstream units without serializing the
    payloads through Python return values.

``permanent``
    Another :class:`.StorageManager` targeting long–term storage. Results saved here
    survive beyond the life of the :class:`.ProtocolDAG` run (for example for
    inspection or for reuse in future extensions).

``stdout`` / ``stderr``
    Optional directories where the engine captures subprocess output triggered
    by the ``ProtocolUnit``.  The directories are removed automatically when the context
    closes.

Keeping these handles bundled together and managed by a context manager lets
``ProtocolUnit`` implementers focus on domain logic while the engine ensures
storage gets flushed and temporary directories are cleaned.


Lifecycle
---------

``Context`` implements a `Python context manager`_.
When the context is exited, ``shared`` and ``permanent`` storage managers flushed tracked files back to their underlying :class:``ExternalStorage``.
Any ``stdout`` or ``stderr`` capture directories are also removed.

.. _Python context manager: https://docs.python.org/3/reference/datamodel.html#context-managers

This means each ``ProtocolUnit``'s ``shared`` and ``permanent`` object are not paths, and should not be treated as such.
Both of these are registries that track if a file should be transferred from its location in ``scratch`` to its final location after completing a unit.

If you want to use some from ``shared`` or ``permanent``, you can use ``ctx.shared.load`` or ``ctx.permanent.load``.
This will allow your unit to fetch those objects from their storage for use.


Using Context inside ProtocolUnits
----------------------------------

Every :meth:`.ProtocolUnit._execute` definition must accept ``ctx: Context`` as its
first argument.  Typical usage looks like the example below.

.. code-block:: python

    from gufe import ProtocolUnit, Context

    class SimulationUnit(ProtocolUnit):

        @staticmethod
        def _execute(ctx: Context, *, setup_result, lambda_window, settings):
            scratch_path = ctx.scratch / f"lambda_{lambda_window}"
            scratch_path.mkdir(exist_ok=True)

            # Read upstream artifacts from ctx.shared
            system_file = ctx.shared.load(setup_result.outputs["system_file"])
            topology_file = ctx.shared.load(setup_result.outputs["topology_file"])


            result_path = ctx.scratch / "some_output.pdb"
            # When you register the filename doesn't matter,
            # just as long as you do it before you return
            result_path_final_location = ctx.permanent.register(result_path)
            # This is an example of running something that you want to save
            simulate(output=result_path)

            # Return only lightweight metadata
            return {
                "lambda_window": lambda_window,
                # We use this because it is already namespaced and can be used between units.
                "result_path": result_path_final_location,
            }

The example above showcases how


Choosing between shared and permanent storage
---------------------------------------------

Both ``ctx.shared`` and ``ctx.permanent`` expose the same :class:`.StorageManager`
API but they serve different audiences:

``ctx.shared``
    Optimized for communication between units in the same DAG execution.  The
    execution backend is free to prune these assets once no downstream unit
    references them.

``ctx.permanent``
    Intended for outputs that should survive beyond the immediate DAG, such as
    user-facing reports or artifacts that will seed future runs.

As a rule of thumb, prefer ``ctx.shared`` unless you have a clear requirement
to keep the data after the ``Protocol`` run concludes.  Small scalar values or
lightweight metadata should still be returned directly from ``_execute`` so
they become part of the ``ProtocolUnitResult`` record.


Interaction with Protocols
--------------------------

``Protocol`` instances do not instantiate ``Context`` directly; they declare ``ProtocolUnit`` objects via ``Protocol._create``.
When an execution backend walks the resulting
``ProtocolDAG`` it constructs a ``Context`` for each unit using the DAG label, unit label, scratch directory, and the configured ``ExternalStorage`` implementations.
The backend might provide different ``ExternalStorage`` implementations (e.g., local filesystem, object store, cluster scratch) depending on where the work runs, but the ``Context`` API seen by ``ProtocolUnit`` authors stays consistent.

Because the execution backend is in charge of creating the contexts, protocol authors can rely on ``ctx`` always being populated with valid storage managers and paths that are safe to write to from distributed workers.


Migrating from legacy Context usage
-----------------------------------

Before ``Context`` was rewritten, it was a simple data class with two ``pathlib.Path`` handles: ``scratch`` and ``shared``.
Existing protocols adopted a variety of implicit conventions around those attributes.
Follow this checklist when migrating old protocols:

1. **Swap file paths to StorageManager APIs.**  Calls like ``ctx.shared /
   "filename"`` should be replaced with the helper methods offered by
   :class:`StorageManager`.  For example:

.. code-block:: python

    path = ctx.shared.scratch_dir / "myfile.dat"
    ctx.shared.register("myfile.dat")


2. **Avoid storing heavy objects in Python outputs.**  Older protocols often
   returned raw ``Path`` objects pointing at scratch files.  Instead, register
   the file with the storage manager and return the storage key (a string) from
   ``_execute`` as shown in the example above.  Downstream units can then call
   :meth:`.StorageManager.load`.

3. **Handle ctx.permanent.**  There was no equivalent in the legacy API.
   Decide which results must persist between DAG executions and write them via
   ``ctx.permanent``.  For migration you can start by mirroring whatever used
   to live in ``ctx.shared`` and refine later.

4. **Expect automatic cleanup.**  Old contexts typically left stdout/stderr
   directories around.  The new context removes these when the unit finishes.
   If your code tried to re-read the capture directories after ``_execute``
   returned, move that logic earlier or rely on the logged data captured in the
   ``ProtocolUnitResult``.

5. **Stop constructing Context manually.**
   Some bespoke execution scripts once instantiated ``Context(scratch=..., shared=...)`` by hand.
   That pattern is obsolete because the constructor now requires ``ExternalStorage`` objects.
   Instead, rely on the execution backend to build contexts.
   For unit tests use the helpers in ``gufe.storage.externalresource`` (e.g., ``MemoryStorage``) to create the necessary storage instances.
