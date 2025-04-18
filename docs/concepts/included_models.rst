Data models included in **gufe**
================================

The core of the **gufe** data model is the :class:`.GufeTokenizable` class,
but **gufe** features more than just this base data structure.

In order to ensure interoperability,
**gufe** also defines classes of objects that represent the core chemistry and alchemistry of a free energy pipeline,
including molecules, chemical systems, and alchemical transformations.
This provides a shared language that tools in the ecosystem use.

Some of these classes are designed to be subclassed, and constitute the *extensible points* of the library.
These include:

1. :ref:`Component <component>`
2. :ref:`Protocol <protocol>`
3. :ref:`AtomMapper <atommapper>`


.. _component:

``Component``
-------------

The :class:`.Component` class is used to represent a portion of a system of molecules,
with a single ``Component`` capable of representing anything from an individual drug-like molecule, an entire protein, or (the concept of) a solvent with ions.

These are often used to define the *components* of a :ref:`chemicalsystem`, which form the nodes of an :ref:`alchemicalnetwork`.
The same ``Component`` may be present within multiple ``ChemicalSystem``\s, such as a :class:`.ProteinComponent` in an ``AlchemicalNetwork`` featuring relative binding transformations between ligands.

As another distinct example: the :class:`.SmallMoleculeComponent` class is used to form the nodes of a :ref:`ligandnetwork`.
This is useful for representing relative transformations between a series of small molecules without invoking the additional complexity of an :ref:`alchemicalnetwork`.

The :class:`.Component` is an *extensible point* of the library,
and is intended to be subclassed to enable new applications.
For details on how to create your own :class:`.Component` classes, see :ref:`howto-component`.


.. _chemicalsystem:

``ChemicalSystem``
------------------

A :class:`.ChemicalSystem` represents a complete system of molecules, and is often composed of multiple :ref:`Components <component>`.


.. _transformation:

``Transformation``
------------------


.. _protocol:

``Protocol``
------------



.. _atommapper:

``AtomMapper``
--------------



.. _ligandnetwork:

``LigandNetwork``
-----------------



.. _alchemicalnetwork:

``AlchemicalNetwork``
---------------------



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

In order to facilitate comparisons of different approaches, **gufe** defines a
hierarchy of simulation settings. This allows certain settings (such as
temperature and pressure) to be consistent across different simulation
tools, while allowing additional custom settings specific to a given tool to
be defined.

Protocols
---------

The actual simulation of a free energy calculation is defined by a **gufe**
:class:`.Protocol`. The :class:`.Protocol` is described as a set of tasks,
each a :class:`.ProtocolUnit`, which may depend on other tasks. As such,
they form a directed acyclic graph.

Gufe does not implement any free energy protocols, but by providing the
abstract API, allows protocol authors to create new simulation protocols
without needing to focus on the details of execution or storage.

Executors
---------

Executors actually run the simulations described by the :class:`.Protocol`.
**gufe** does not define an executor API, although it includes the very simple
serial executor in :func:`.execute_DAG`.

The responsibilities of an executor include running the tasks (units) for a
:class:`.Protocol` and managing storage of output. Gufe contains some tools
to facilitate that, particularly around storage, but it is up to the
executor to determine how to/whether to use those.
