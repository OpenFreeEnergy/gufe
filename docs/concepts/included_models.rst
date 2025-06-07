Data models included in **gufe**
================================

The core of the **gufe** data model is the :class:`.GufeTokenizable` class,
but **gufe** features more than just this base data structure.

In order to ensure interoperability,
**gufe** also defines classes of objects that represent the core chemistry and alchemistry of a free energy pipeline,
including molecules, chemical systems, and alchemical transformations.
This provides a shared language that tools in the ecosystem use.


.. note::

    Some of these classes are designed to be subclassed, and constitute the *extensible points* of the library.
    These include the following; see the **How-To Guide** for more information on how to extend from each:
    
    1. :ref:`Component <component>` : :ref:`howto-component`
    2. :ref:`Protocol <protocol>` : :ref:`howto-protocol`
    3. :ref:`ComponentMapping <componentmapping>` : :ref:`howto-componentmapping`


.. _component:

``Component``
-------------

The :class:`.Component` class represents a portion of a system of molecules,
with a single ``Component`` capable of representing anything from an individual drug-like molecule, an entire protein, or (even the concept of) a solvent with ions.

These are often used to define the *components* of a :ref:`chemicalsystem`, which form the nodes of an :ref:`alchemicalnetwork`.
The same ``Component`` may be present within multiple ``ChemicalSystem``\s, such as a :class:`.ProteinComponent` in an ``AlchemicalNetwork`` featuring relative binding transformations between ligands.

As another distinct example: the :class:`.SmallMoleculeComponent` class is used to form the nodes of a :ref:`ligandnetwork`.
This is useful for representing relative transformations between a series of small molecules without invoking the additional complexity of an :ref:`alchemicalnetwork`.

.. note::
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

When used as inputs to a ``Transformation``, ``ChemicalSystem``\s represent the set of ``Component``\s for which a free energy difference will be estimated.
Alchemical methods performing free energy perturbation (FEP) between the two ``ChemicalSystem``\s of a ``Transformation`` will simulate these ``Component``\s using some sampling approach, obtaining enough information to derive a free energy difference estimate.


.. _transformation:

``Transformation``
------------------

A :class:`.Transformation` represents an alchemical transformation between two :ref:`ChemicalSystems <chemicalsystem>`.

``Transformation`` objects are often used as the edges of an :ref:`alchemicalnetwork`.
In addition to referencing the ``ChemicalSystem``\s it spans,
a ``Transformation`` also includes the :ref:`protocol` used to actually perform the alchemical transformation,
as well as an :ref:`componentmapping` specifying what portions of the :ref:`Components <component>` are being transformed across the ``ChemicalSystem``\s.

A ``Transformation`` functions as a container for all the information needed to obtain an estimate of the free energy difference between its two ``ChemicalSystem``\s.


.. _nontransformation:

``NonTransformation``
---------------------

A :class:`.NonTransformation` represents non-alchemical sampling of a single :ref:`ChemicalSystem <chemicalsystem>`.

In the context of an :ref:`alchemicalnetwork`, a ``NonTransformation`` is effectively a self-loop, featuring the same ``ChemicalSystem`` at either end..
Similar to a :ref:`Transformation <transformation>`, it features a :ref:`protocol` used to perform sampling on its ``ChemicalSystem``, but does not feature a :ref:`componentmapping` since none is required for this.
An example of a ``Protocol`` that would be appropriate for a ``NonTransformation`` is one that performs equilibrium molecular dynamics of the ``ChemicalSystem``.

A ``NonTransformation`` cannot be used to obtain a free energy difference estimate, since by definition transforming the ``ChemicalSystem`` to itself should give exactly ``0``.


.. _protocol:

``Protocol``
------------

A :class:`.Protocol` represents the specific sampling approach used to transform one :ref:`ChemicalSystem <chemicalsystem>` into another (as in a :ref:`Transformation <transformation>`), or to simply sample a single :ref:`ChemicalSystem <chemicalsystem>` (as in a :ref:`NonTransformation <nontransformation>`).

``Protocol`` objects are often used as part of a ``Transformation``, although they can be used on their own alongside ``ChemicalSystem``\s and ``ComponentMapping``\s (when needed) to obtain free energy difference estimates.
Individual ``Protocol`` subclasses obtain these estimates in a wide variety of ways, with varying domains of applicability and effectiveness.

.. note::
    The :class:`.Protocol` is an *extensible point* of the library,
    and is intended to be subclassed to enable new applications.
    For details on how to create your own :class:`.Protocol` classes, see :ref:`howto-protocol`.

.. _protocoldag:

``ProtocolDAG``
^^^^^^^^^^^^^^^

A :class:`.ProtocolDAG` is an executable artifact that performs a :ref:`Protocol <protocol>`.

A ``ProtocolDAG`` is created via :meth:`.Protocol.create` in combination with :ref:`ChemicalSystem(s) <chemicalsystem>` and a :ref:`ComponentMapping <componentmapping>` (when needed). 
It is a `directed acyclic graph <https://en.wikipedia.org/wiki/Directed_acyclic_graph>`_ (DAG) of :ref:`ProtocolUnits <protocolunit>` and their dependency relationships.
The ``ProtocolUnit``\s of this ``ProtocolDAG`` can be executed in dependency-order to yield information needed for a free energy difference estimate.

``ProtocolDAG``\s are generally only handled directly by ecosystem tools that perform :ref:`Transformation <transformation>` execution.


.. _protocolunit:

``ProtocolUnit``
^^^^^^^^^^^^^^^

A :class:`.ProtocolUnit` is the unit of execution of a :ref:`ProtocolDAG <protocoldag>`, functioning as a node with dependency relationships in the DAG.

A ``ProtocolUnit`` 


.. note::
    The :class:`.Protocol` is an *extensible point* of the library,
    and is intended to be subclassed to enable new applications.
    For details on how to create your own :class:`.Protocol` classes, see :ref:`howto-protocol`.



.. _protocolunitresult:

``ProtocolUnitResult``
^^^^^^^^^^^^^^^^^^^^^^

A :class:`.ProtocolUnitResult` is created by a :ref


.. _protocolunitfailure:

``ProtocolUnitFailure``
^^^^^^^^^^^^^^^^^^^^^^^

A :class:`.ProtocolUnitFailure` is created by a :ref


.. _protocoldagresult:

``ProtocolDAGResult``
^^^^^^^^^^^^^^^^^^^^^

A :class:`.ProtocolDAGResult` retains the results from executing a :ref:`ProtocolDAG <protocoldag>`.

A ``ProtocolDAGResult`` contains the same information as a ``ProtocolDAG`` (including ``ProtocolUnit``\s and their dependency relationships), while also featuring the set of :ref:`ProtocolUnitResults <protocolunitresult>` that resulted from each.
``ProtocolDAGResult``\s form the smallest 


.. _protocolresult:

``ProtocolResult``
^^^^^^^^^^^^^^^^^^

A :class:`.ProtocolResult` aggregates the results from one or more :ref:`ProtocolDAGResults <protocoldagresult>` to yield a free energy difference estimate.

``ProtocolResult`` are created from :meth:`.Protocol.gather`, and feature the ``Protocol``-specific methods necessary to obtain actual free energy difference estimates from a set of ``ProtocolDAGResult``\s, namely:

* :meth:`.ProtocolResult.get_estimate`
* :meth:`.ProtocolResult.get_uncertainty`


.. note::
    The :class:`.ProtocolResult` is an *extensible point* of the library alongside :class:`Protocol`,
    and is intended to be subclassed to enable new applications.
    For details on how to create your own :class:`.ProtocolResult` classes, see :ref:`howto-protocol`.



.. _componentmapping:

``ComponentMapping``
--------------------

.. note::
    The :class:`.ComponentMapping` is an *extensible point* of the library,
    and is intended to be subclassed to enable new applications.
    For details on how to create your own :class:`.ComponentMapping` classes, see :ref:`howto-componentmapping`.


.. _ligandnetwork:

``LigandNetwork``
-----------------



.. _alchemicalnetwork:

``AlchemicalNetwork``
---------------------



