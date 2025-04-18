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

A :class:`.ChemicalSystem` represents a complete system of molecules,
and is often composed of multiple :ref:`Components <component>`.

These are most often used as nodes of an :ref:`alchemicalnetwork`, with pairs of ``ChemicalSystem``\s connected by :ref:`Transformations <transformation>`.
Because a ``ChemicalSystem`` functions as a kind of container of :ref:`Components <component>`, more than one ``ChemicalSystem`` can feature the same ``Component``\s.
This allows even very large ``AlchemicalNetwork``\s to be relatively small in memory, as only a few large ``Component``\s like :class:`.ProteinComponent`\s may be shared among hundreds of ``ChemicalSystem``\s.


.. _transformation:

``Transformation``
------------------

A :class:`.Transformration` represents an alchemical transformation between two :ref:`ChemicalSystems <chemicalsystem>`.




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



