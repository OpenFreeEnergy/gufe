
Understanding GufeTokenizables
==============================

Most objects in gufe are subclasses of :class:`.GufeTokenizable`.
This base class enforces common behavior and sets requirements necessary
to guarantee performance and reproducibility between all downstream packages that use gufe.

For example, when we create a ``SmallMoleculeComponent`` representing benzene, that object is also a ``GufeTokenizable``:

.. code:: python

    >>> import gufe
    >>> benzene = gufe.SmallMoleculeComponent.from_sdf_file("benzene.sdf")
    >>> type(Benzene)
    gufe.components.smallmoleculecomponent.SmallMoleculeComponent

    >>> from gufe.tokenization import GufeTokenizable
    >>> isinstance(benzene, GufeTokenizable)
    True

By definition, a ``GufeTokenizable`` must be:

1. :ref:`immutable <immutability>`
2. :ref:`hashable <gufe_keys>`
3. :ref:`serializable <serialization>`

.. _immutability:

1. Immutability of ``GufeTokenizables``
---------------------------------------

One important restriction on :class:`.GufeTokenizable` subclasses is that they must be immutable,
meaning that none of its attributes change after initialization.
In other words, all attributes should be set when you create an object, and never changed after that.
If your object is immutable, then it is suitable to be a GufeTokenizable.

For example, once the benzene molecule from above is loaded, its attributes (such as ``name``) are immutable:

.. code:: python

    >>> benzene.name
    'benzene'
    >>> benzene.name = 'benzene_1'
    AttributeError

.. TODO: note that no error is raised if we try to mutate the dict object, e.g. ``benzene.to_dict()['atoms'] = 1``?

Immutability is critical to gufe's design, because it means that gufe can generate a deterministic unique identifier (the gufe key)
based on the ``GufeTokenizable``'s properties.


.. TODO: talk about `copy_with_replacements`?

.. TODO: how to actually implement a mutable attribute? isn't this enforced, or does this just mean using the unfreeze functionality?

.. There is a special case of mutability that is also allowed, which is if the
.. object is functionally immutable.  As an example, consider a flag to turn on
.. or off usage of a cache of input-output pairs for some deterministic method.
.. If the cache is turned on, you first try to return the value from it, and
.. only perform the calculation if the inputs don't have a cached output
.. associated. In this case, the flag is mutable, but this has no effect on the
.. results. Indeed, the cache itself may be implemented as a mutable attribute
.. of the object, but again, this would not change the results that are
.. returned. It would also be recommended that an attribute like a cache, which
.. is only used internally, should be marked private with a leading underscore.
.. On the other hand, a flag that changes code path in a way that might
.. change the results of any operation would mean that the object cannot be a
.. :class:`.GufeTokenizable`.

.. _gufe_keys:

2. Hashing ``GufeTokenizables``: the gufe key
---------------------------------------------

Because gufe objects are immutable, each object has a unique identifier, which we call its ``key``.
The ``key`` is a string, typically in the format ``{CLASS_NAME}-{HEXADECIMAL_LABEL}``.

For our benzene ``SmallMoleculeComponent``, the key is ``'SmallMoleculeComponent-ec3c7a92771f8872dab1a9fc4911c795``:

.. code:: python

    >>> benzene.key
    'SmallMoleculeComponent-ec3c7a92771f8872dab1a9fc4911c795'

For most objects, the hexadecimal label is generated based on the contents of the class -- in
particular, it is based on contents of the ``_to_dict()`` dictionary, filtered
to remove anything that matches the ``_defaults()`` dictionary.

This gives the gufe key the following important properties:

* A key is based on a **cryptographic hash**, so it is extremely unlikely
  that two objects that are functionally different will have the same key.
* Key creation is **deterministic**, so that it is preserved across different creation times,
  including across different hardware, across different Python sessions,
  and even within the same Python session.
* A key is preserved across minor versions of the code, since it is dependent on non-default attributes and
   we follow `SemVer <https://semver.org>`_.

..  QUESTION: is this still true, or have we changed keys across minor versions?

These properties, in particular the stability across Python sessions,  make the gufe key a stable identifier for the object.
This stability means that they can be used for store-by-reference, and therefore deduplicated to optimize memory and performance.

Deduplication of GufeTokenizables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are two types of deduplication of GufeTokenizables.
Objects are deduplicated in memory because gufe keeps a registry of all instantiated GufeTokenizables.
Objects can be deduplicated on storage to disk because we store by reference to the gufe key. 

.. _gufe-memory-deduplication:

Deduplication in memory (flyweight pattern)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Memory deduplication means that only one object with a given gufe ``key``
will exist in any single Python session. 
We ensure this by maintaining a registry of all GufeTokenizables that gets updated any time a
GufeTokenizable is created. (The registry is a mapping to weak references, which
allows Python's garbage collection to clean up GufeTokenizables that are no
longer needed.) This is essentially an implementation of the `flyweight
pattern <https://en.wikipedia.org/wiki/Flyweight_pattern>`_.

This memory deduplication is ensured by the ``GufeTokenizable.from_dict``,
which is typically used in deserialization. It will always use the first
object in memory with that ``key``. This can lead to some unexpected
behavior; for example, using the ``Foo`` class defined above:

.. code::

    # here Foo is a GufeTokenizable:
    >>> a = Foo(0)
    >>> b = Foo(0)
    >>> a is b
    False
    >>> c = Foo.from_dict(a.to_dict())
    >>> c is a  # surprise!
    True
    >>> d = Foo.from_dict(b.to_dict())
    >>> d is b
    False
    >>> d is a  # this is because `a` has the spot in the registry
    True


Deduplication on disk
~~~~~~~~~~~~~~~~~~~~~

Deduplication in disk storage is fundamentally the responsibility of the
specific storage system, which falls outside the scope of ``gufe``.
However, ``gufe`` provides some tools to facilitate implementation of a storage
system.

The main idea is to use the ``key`` to ensure uniqueness, and to use it as a label for the object's serialized representation.
Additionally, the ``key``, which is simply a string, can be used as a stand-in for the object.
When an outer GufeTokenizable contains an inner GufeTokenizable, the outer can store the key in place of the inner object.
That is, we can store by reference to the key.

To convert a GufeTokenizable ``obj`` into a dictionary that references inner
GufeTokenizables by key, use ``obj.to_keyed_dict()``. That method replaces
each GufeTokenizable by a dict with a single key, ``':gufe-key:'``, mapping
to the key of the object. Of course, you'll also need to do the same for all
inner GufeTokenizables; to get a list of all of them, use
:func:`.get_all_gufe_objs` on the outermost ``obj``.

.. TODO: add a tutorial for this


.. _serialization:

3. Serialization (and deserialized representations)
---------------------------------------------------

Any GufeTokenizable can represented in the following ways:

.. find nice simple but nested test data to demo this


a) dictionary
^^^^^^^^^^^^^

The ``to_dict()`` method is the most explicit way to represent a GufeTokenizable. 
This method recursively unpacks any inner GufeTokenizables that an
outer GufeTokenizable contains to their full dict representation.
Although this method is best way to see all information stored in a GufeTokenizable,
it is also the least space-efficient.

.. TODO: show this method
.. TODO: diagram

b) shallow dictionary
^^^^^^^^^^^^^^^^^^^^^

The ``to_shallow_dict()`` method is similar to ``to_dict()`` in that it unpacks a tokenizable into a ``dict`` format,
but a shallow dict is *not recursive* and only unpacks the top level of the GufeTokenizable. Anything nested deeper is represented by
the inner objects' GufeTokenizable.

.. TODO: show this method
.. TODO: diagram


This method is most useful for iterating through the hierarchy of a GufeTokenizable one layer at a time.


c) keyed dictionary
^^^^^^^^^^^^^^^^^^^

The ``to_keyed_dict()`` method is similar to ``to_shallow_dict`` in that it only unpacks the first layer of a GufeTokenizable.
However, a keyed dict represents the next layer as its gufe key, e.g. ``{':gufe-key:': 'ChemicalSystem-96f686efdc070e01b74888cbb830f720'},``
  
A keyed dict is the most compact representation of a GufeTokenizable and can be useful for understanding its contents,
but it does not have the complete representation for reconstruction or sending information (for this, see the next section, :ref:`keyed chain <keyed_chain>`)

.. TODO: show this method
.. TODO: diagram

.. _keyed_chain:

d) keyed chain
^^^^^^^^^^^^^^

The ``keyed_chain()`` method is a powerful representation of a GufeTokenizable that enables efficient reconstruction of an object without duplication.
It uses ``keyed_dict`` to unpack a GufeTokenizable from the bottom (innermost) layer up, effectively constructing a DAG
(`directed acyclic graph <https://en.wikipedia.org/wiki/Directed_acyclic_graph>`_) where re-used GufeTokenizables are deduplicated.

.. TODO: maybe show output, maybe abbreviated?
.. TODO: diagram (especially this one!!)

.. NOTE::
  See :doc:`../how-tos/serialization` for details on how to implement serialization of your own GufeTokenizables.