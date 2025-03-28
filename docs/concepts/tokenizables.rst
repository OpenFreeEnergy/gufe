
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

1. Immutability of GufeTokenizables
-----------------------------------

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

.. todo: note that no error is raised if we try to mutate the dict object, e.g. ``benzene.to_dict()['atoms'] = 1``?

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

2. Hashing GufeTokenizables: the gufe key
-----------------------------------------

.. TODO: code snippet showing key

Because gufe objects are immutable, each object has a unique identifier, which we call its ``key``.
The ``key`` is a string, typically in the format ``{CLASS_NAME}-{HEXADECIMAL_LABEL}``.

For our benzene ``SmallMoleculeComponent``, the key is ``'SmallMoleculeComponent-ec3c7a92771f8872dab1a9fc4911c795``:

.. code:: python

    >>> benzene.key
    'SmallMoleculeComponent-ec3c7a92771f8872dab1a9fc4911c795'

.. code snippet showing object instantiation and key
.. maybe also show how we can manually create the key and it's just based on the filtered dict?

For most objects, the hexadecimal label is generated based on the contents of the class -- in
particular, it is based on contents of the ``_to_dict()`` dictionary, filtered
to remove anything that matches the ``_defaults()`` dictionary.

This gives the gufe key the following important properties:

* A key is based on a **cryptographic hash**, so it is extremely unlikely
  that two objects that are functionally different will have the same key.
* Key creation is **deterministic**, so that it is preserved across different creation times,
  including across different hardware, across different Python sessions,
  and even within the same Python session.
* A key is based on the non-default attributes, it is preserved across minor versions of the code
  (since we follow `SemVer <https://semver.org>`_).

These properties, in particular the stability across Python sessions,  make the gufe key a stable identifier for the object.
This stability means that they can be used for store-by-reference, and therefore deduplicated to optimize memory and performance.

Deduplication of GufeTokenizables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are two types of deduplication of GufeTokenizables. Objects can be
deduplicated on storage to disk because we store by reference to the gufe
key. Additionally, objects are deduplicated in memory because we keep a
registry of all instantiated GufeTokenizables.

.. _gufe-memory-deduplication:

Deduplication in memory (flyweight pattern)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Memory deduplication means that only one object with a given gufe ``key``
will exist in any single Python session. We ensure this by maintaining a
registry of all GufeTokenizables that gets updated any time a
GufeTokenizable is created. (This is a mapping to weak references, which
allows Python's garbage collection to clean up GufeTokenizables that are no
longer needed.) This is essentially an implementation of the `flyweight
pattern <https://en.wikipedia.org/wiki/Flyweight_pattern>`_.

This memory deduplication is ensured by the ``GufeTokenizable.from_dict``,
which is typically used in deserialization. It will always use the first
object in memory with that ``key``. This can lead to some unexpected
behavior; for example, using the ``Foo`` class defined above:

.. code::

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
specific storage system, which falls outside the scope of ``gufe``. However,
``gufe`` provides some tools to facilitate implementation of a storage
system.

The main idea is again to use the ``key`` to ensure uniqueness, and to use
it as a label for the object's serialized representation.  Additionally, the
``key``, as a simple string, can be used as a stand-in for the object, so
when an outer GufeTokenizable contains an inner GufeTokenizable, the
outer can store the key in place of the inner object.  That is, we can store
by reference to the key.

To convert a GufeTokenizable ``obj`` into a dictionary that references inner
GufeTokenizables by key, use ``obj.to_keyed_dict()``. That method replaces
each GufeTokenizable by a dict with a single key, ``':gufe-key:'``, mapping
to the key of the object. Of course, you'll also need to do the same for all
inner GufeTokenizables; to get a list of all of them, use
:func:`.get_all_gufe_objs` on the outermost ``obj``.


.. _serialization:

3. Serialization
----------------

Any GufeTokenizable can represented in the following ways:

.. find nice simple but nested test data to demo this


1. dict

  - this is the most explicit way to represent a GufeTokenizable
  - unpacks all levels to their dict representation
  - complete but not space efficient

2. shallow_dict

  - only one level is 'unpacked', anything deeper is stored by GufeTokenizable
  - (QUESTION: where/how does this happen?, `to_shallow_dict` is just calling to_dict?)
  - most useful for iterating through a gufe tokenizable layer-by-layer

3. keyed_dict

  - similar to shallow_dict, only one level is unpacked, anything deeper is represented as
    {':gufe-key:': 'ChemicalSystem-96f686efdc070e01b74888cbb830f720'},
  - most compact representation of the object, but does not have the complete representation for serialization, sending information

4. keyed_chain

  - uses keyed_dict to create a DAG for efficient reconstruction without duplication
  - explain with a diagram here?

See :doc:`../how-tos/serialization` for details on how to implement serialization of your own GufeTokenizables.