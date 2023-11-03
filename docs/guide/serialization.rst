Serializing GUFE objects
========================

When you create objects the follow the GUFE APIs, you will want them to be
serializable; that is, you'll want to be able to store everything to disk.
In this section, we'll discuss how GUFE handles serialization, and what a
developer needs to do to serialize an object.

From a user perspective, this information is mainly relevant for
understanding what you can expect from storage systems built for GUFE
objects. From a developer perspective, this will tell you what you need to
do to be compatible with (and get the most from) the GUFE serialization
mechanisms.

In practice: What a developer needs to do
-----------------------------------------

For most use cases, an object can be made serializable by gufe if it (1)
inherits from :class:`.GufeTokenizable` and (2) implements several abstract
private methods related to serialization: ``_to_dict``, ``_from_dict``, and
``_defaults``:

* ``_to_dict`` returns a dictionary mapping string names for each attribute
  to preserve to a serializable representation of that attribute. A
  representation is considered "serializable" if it is (1) a
  :class:`.GufeTokenizable`, (2) a primitive type that can be handled by
  Python's built-in JSON, or (3) a type for which a custom
  :class:`.JSONCodec` has been registered (see :ref:`customjson`).
* ``_from_dict`` is a ``classmethod`` that takes the dictionary returned by
  ``_to_dict`` and reconstructs the object.
* ``_defaults`` returns whatever ``_to_dict`` would return for any default
  values; i.e., it should be a dictionary where the keys are the the keys of
  the ``_to_dict`` dictionary that have default values, and the values are
  the associated defaults. Importantly, this may differ from the default
  value given at initialization -- for example, it is common practice to use
  ``None`` as a stand-in for a mutable container (e.g., an empty list). In
  that case, the value given by ``_defaults`` should be the container, not
  ``None``, since that is what will be in the ``_to_dict``.

As a very simple example, consider the following complete implementation:

.. code::

    class Foo(GufeTokenizable):
        def __init__(self, bar, baz=None):
            self.bar = bar
            if baz is None:
                baz = []
            self.baz = baz

        def _to_dict(self):
            return {'bar': self.bar, 'baz': self.baz}

        @classmethod
        def _from_dict(cls, dct):
            return cls(**dct)

        def _defaults(self):
            return {'baz': []}

Testing your GufeTokenizable
----------------------------

We provide a convenience mix-in for testing that a
:class:`.GufeTokenizable` will be correctly stored and reloaded. Here's an
example of how to use it for the ``Foo`` class above:

.. code::

    import pytest
    from gufe.tests.test_tokenization import GufeTokenizationTestsMixin

    class TestFoo(GufeTokenizationTestsMixin):
        cls = Foo
        key = "Foo-8e7fd2803d73ef0cd9faceef751acfca"

        @pytest.fixture
        def instance(self):
            return Foo(5, baz=['qux', 'quux', 'quuux'])

You'll need to get the ``key`` from creating the ``instance`` once and
recording ``str(obj.key)``-- since the key is stable, it shouldn't change
(and, in fact, that is tested). This will add several tests to your test
suite to ensure that you can save and reload your :class:`.GufeTokenizable`.

.. TODO: add a section here about the various to_*_dict from_*_dict forms:
   shallow_dict; (deep) dict; keyed_dict

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

.. _customjson:

Supporting custom types in serialization
----------------------------------------

Occasionally, you may need to add support for a custom type that is not
already known by GUFE. For example, there may be objects that are included
in your GUFE object that come from another library, and are not GUFE
tokenizables.

When this happens, you will need to create and register a custom JSON
codec. To do this, at a minimum you need to define a ``to_dict`` and
``from_dict`` function for your type. You can then give the class of your
type and those functions to a :class:`.JSONCodec`.  As an example, consider
the code for our codec for ``pathlib.Path`` objects:

.. code::

    PATH_CODEC = JSONCodec(
        cls=pathlib.Path,
        to_dict=lambda p: {'path': str(p)},
        from_dict=lambda dct: pathlib.Path(dct['path'])
    )

Here the ``to_dict`` and ``from_dict`` are lambdas. The ideas of these are
the same as the ``GufeTokenizable._to_dict`` and
``GufeTokenizable._from_dict`` described above.

In this default setup, the codec recognizes your object type by doing an
``isinstance`` check on the ``cls`` you gave it. It updates the dictionary
that comes from your ``to_dict`` with the class and module of the object,
as well as a key to mark this dictionary as coming from this codec. The
key-value pairs that it adds make it so that the codec can recognize the
dictionary when it is deserializes (decodes) the data.

Details of how the object is recognized for encoding or how the dictionary
is recognized for decoding can be changed by passing functions to the
``is_my_obj`` or ``is_my_dict`` parameters of :class:`.JSONCodec`.

.. warning::
    The Custom encoders & decoders only override the default JSON
    encoder/decoder if they are not able to natively handle the object.
    This leads to some odd / lossy behaviour for some objects such
    as ``np.float64`` which is natively converted to a ``float`` type
    by the default encoder, whilst other numpy generic types are
    appropriately roundtripped.

On the use of NumPy generic types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Due to their inconsistent behaviour in how they are handled by the default
JSON encoder/decoder routines (see warning above), it is our suggestion
that Python types should be used preferrentially instead of NumPy generic
types. For example if one would be looking to store a single float value,
a ``float`` would be prefered to a ``np.float32`` or ``np.float64``.

Please note that this only applied to generic types being used **outside of
numpy arrays**. NumPy arrays are, as far as we know, always handled
in a consistent manner.

Dumping arbitrary objects to JSON
---------------------------------

Any :class:`.GufeTokenizable` can be dumped to a JSON file using the custom
JSON handlers. Given a :class:`.GufeTokenizable` called ``obj`` and a
path-like called ``filename``, you can dump to JSON with this recipe:

.. code::

    import json
    from gufe.tokenization import JSON_HANDLER
    with open(filename, mode='w') as f:
        json.dump(obj.to_dict(), f, cls=JSON_HANDLER.encoder)

Similarly, you can reload the object with:

.. code::

    import json
    from gufe.tokenization import JSON_HANDLER
    with open(filename, mode='r') as f:
        obj = json.load(f, cls=JSON_HANDLER.decoder)

Note that these objects are not space-efficient: that is, if you have
the same object in memory referenced by multiple objects (e.g., an identical
``ProteinComponent`` in more than one ``ChemicalSystem``), then you will
save multiple copies of its JSON representation.

On reloading, tools that use the recommended ``from_dict`` method will undo
do this duplication; see :ref:`gufe-memory-deduplication` for details.

.. Using JSON codecs outside of JSON
.. ---------------------------------

.. In a custom recursive storage scheme, 

.. TODO: DWHS wants to write something here that describes how to use the
   codecs in your own non-JSON storage scheme. But this is complicated
   enough that it will take significant time

When your object has recursive references
-----------------------------------------

In some cases, your object may have recursive references to other objects.
For example, you may have objects ``parent`` and ``child``, where
``parent.children`` includes ``child`` and ``child.parent`` is ``parent``.
This means that the object dependency graph is not a directed acyclic graph.
The best solution here is to avoid this design pattern whenever possible.
Importantly, the ``child`` object in this case cannot be immutable, unless
it is only created as a part of the creation of ``parent``.

However, this could be made functionally immutable by (1) requiring that a
valid ``child`` set its ``parent`` attribute exactly once; (2) not storing
the ``child.parent`` attribute, and instead ensuring it gets set by the
``parent`` (e.g., in ``__init__``). Note that this also assumes that all
``children`` are provided to ``parent`` on initialization, otherwise
``parent`` is also mutable.

Another approach here would be to make it so that ``child`` was not a
:class:`.GufeTokenizable`, and instead explicitly handle its serialization
in a more complicated ``parent._to_dict`` method, with similarly more
complicated ``parent._from_dict``. This approach is particularly useful if
the ``child`` object isn't too complicated, and if the ``child`` is unlikely
to be reused outside the context of the ``parent``.

Understanding the theory: The GUFE key
--------------------------------------

One of the important concepts is that every GUFE object has a unique
identifier, which we call its ``key``. The ``key`` is a string, typically
in the format ``{CLASS_NAME}-{HEXADECIMAL_LABEL}``, e.g.,
``ProteinComponent-7338abda590510f1dae764e068a65fdc``. For most objects, the
hexadecimal label is generated based on the contents of the class -- in
particular, it is based on contents of the ``_to_dict`` dictionary, filtered
to remove anything that matches the ``_defaults`` dictionary.

This gives the GUFE key a number of important properties:

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
the GUFE key a stable identifier for the object, which means that they can
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

Memory deduplication means that only one object with a given GUFE ``key``
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
