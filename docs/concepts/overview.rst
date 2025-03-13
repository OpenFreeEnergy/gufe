Overview
========

Gufe exists to define interoperable APIs for free energy calculations, which
can be used by an ecosystem of Python packages to develop tools that focus
on performing specific aspects of the free energy pipeline, while
benefitting from other tools for other aspects.

In order to do this, gufe distinguishes several aspects of the process of
calculating free energies, and define APIs for each of them. Gufe provides
the underlying infrastructure; the real science is done in packages that
implement parts of the gufe APIs.

Core data models
----------------

In order to ensure interoperability, gufe defines objects that represent the
core chemistry and alchemistry of a free energy pipeline, including
molecules, chemical systems, and alchemical transformations. This provides a
shared language that different tools use.

Ligand network setup
--------------------

Gufe defines a basic API for the common case of performing alchemical
transformations between small molecules, either for relative binding free
energies of relative hydration free energies. This handles how mappings
between different molecules are defined for alchemical transformations,
by defining both the :class:`.LigandAtomMapping` object that contains the
details of a specific mapping, and the :class:`.AtomMapper` abstract API for
an object that creates the mappings.

Simulation settings
-------------------

In order to facilitate comparisons of different approaches, gufe defines a
hierarchy of simulation settings. This allows certain settings (such as
temperature and pressure) to be consistent across different simulation
tools, while allowing additional custom settings specific to a given tool to
be defined.

Protocols
---------

The actual simulation of a free energy calculation is defined by a gufe
:class:`.Protocol`. The :class:`.Protocol` is described as a set of tasks,
each a :class:`.ProtocolUnit`, which may depend on other tasks. As such,
they form a directed acyclic graph.

Gufe does not implement any free energy protocols, but by providing the
abstract API, allows protocol authors to create new simulation protocols
without needing to focus on the details of execution or storage.

Executors
---------

Executors actually run the simulations described by the :class:`.Protocol`.
gufe does not define an executor API, although it includes the very simple
serial executor in :func:`.execute_DAG`.

The responsibilities of an executor include running the tasks (units) for a
:class:`.Protocol` and managing storage of output. Gufe contains some tools
to facilitate that, particularly around storage, but it is up to the
executor to determine how to/whether to use those.

Strategies
----------

Strategies have yet to be implemented, but the gufe design leaves a place
for an object that, at the scale of an alchemical network, can dynamically
decide where to focus more simulation effort based on the results that have
been received so far. This will be useful for adaptive approaches to
sampling a network.

Core gufe infrastructure
------------------------

Behind the scenes, gufe implements a number of details common to all of its
objects. Nearly all objects in gufe are subclasses of
:class:`.GufeTokenizable`, which sets a few requirements and behaviors on
these objects, including that they are immutable, serializable, and
hashable. These provide important guarantees to downstream packages that
facilitate repeatability and reproducibility.
