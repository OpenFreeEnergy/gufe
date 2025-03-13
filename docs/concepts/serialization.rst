
GufeTokenizables must be functionally immutable
-----------------------------------------------

One important restriction on :class:`.GufeTokenizable` subclasses is that
they must be functionally immutable. That is, they must have no mutable
attributes that change their functionality. In most cases, this means that
the object must be strictly immutable.

When an object is immutable, that means that none of its attributes change
after initialization. So all attributes should be set when you create an
object, and never changed after that. If your object is immutable, then it
is suitable to be a :class:`.GufeTokenizable`.

There is a special case of mutability that is also allowed, which is if the
object is functionally immutable.  As an example, consider a flag to turn on
or off usage of a cache of input-output pairs for some deterministic method.
If the cache is turned on, you first try to return the value from it, and
only perform the calculation if the inputs don't have a cached output
associated. In this case, the flag is mutable, but this has no effect on the
results. Indeed, the cache itself may be implemented as a mutable attribute
of the object, but again, this would not change the results that are
returned. It would also be recommended that an attribute like a cache, which
is only used internally, should be marked private with a leading underscore.

On the other hand, a flag that changes code path in a way that might
change the results of any operation would mean that the object cannot be a
:class:`.GufeTokenizable`.




.. Using JSON codecs outside of JSON
.. ---------------------------------

.. In a custom recursive storage scheme,

.. TODO: DWHS wants to write something here that describes how to use the
   codecs in your own non-JSON storage scheme. But this is complicated
   enough that it will take significant time



Understanding the theory: The gufe key
--------------------------------------

One of the important concepts is that every gufe object has a unique
identifier, which we call its ``key``. The ``key`` is a string, typically
in the format ``{CLASS_NAME}-{HEXADECIMAL_LABEL}``, e.g.,
``ProteinComponent-7338abda590510f1dae764e068a65fdc``. For most objects, the
hexadecimal label is generated based on the contents of the class -- in
particular, it is based on contents of the ``_to_dict`` dictionary, filtered
to remove anything that matches the ``_defaults`` dictionary.

This gives the gufe key a number of important properties:

* Because the key is based on a cryptographic hash, it is extremely unlikely
  that two objects that are functionally different will have the same key.
* Because the key is created deterministically, it is preserved across
  different creation times, including across different hardware, across
  different Python sessions, and even within the same Python session.
* Because the key is based on the non-default attributes, it is preserved
  across minor versions of the code (since we follow `SemVer
  <https://semver.org>`_). The developer-provided ``_defaults`` method is
  used here.

These properties, in particular the stability across Python sessions,  make
the gufe key a stable identifier for the object, which means that they can
be used for store-by-reference.

Deduplication of GufeTokenizables
---------------------------------

There are two types of deduplication of GufeTokenizables. Objects can be
deduplicated on storage to disk because we store by reference to the gufe
key. Additionally, objects are deduplicated in memory because we keep a
registry of all instantiated GufeTokenizables.

.. _gufe-memory-deduplication:

Deduplication in memory
~~~~~~~~~~~~~~~~~~~~~~~

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
