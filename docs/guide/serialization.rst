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

GufeTokenizables must be functionally immutable
-----------------------------------------------

One important restriction on :class:`.GufeTokenizable` subclasses is that
they must be functionally immutable. That is, they must have no mutable
attributes that change their functionality. In most cases, this means that
the object must be strictly immutable.

As an example of a mutable attribute that does not change functionality,
consider a flag to turn on or off usage of a cache of input-output pairs for
some deterministic method. If the cache is turned on, you first try to
return the value from it, and only perform the calculation if the inputs
don't have a cached output associated. In this case, the flag is mutable,
but this has no effect on the results.  Indeed, the cache itself may be
implemented as a mutable attribute of the object, but again, this would not
change the results that are returned.

On the other hand, a flag that changes code path in a way that might
change the results of any operation would not be allowed.

.. _customjson:

Supporting custom types in serialization
----------------------------------------

Occasionally, you may need to add support for a custom type that is not
already known by GUFE. For example, there may be objects that are included
in your GUFE object that come from another library, and are not GUFE
tokenizables.

When this happens, you will need to create and register a custom JSON
codec. To do this, at a minimum you need to define a ``to_dict`` and
``from_dict`` method for your type. Additionally, the  ???

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
the same object in memory stored at multiple locations (e.g., an identical
``SolventComponent`` in more than one ``ChemicalSystem``), then you will
save multiple copies of its JSON representation.

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
in the format ``{CLASS_NAME}-{HEXADECIMAL_LABEL}``. For most objects, the
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
  across minor versions of the code.

These properties make the GUFE key a stable identifier for the object, which
means that they can be used for store-by-reference.  When one
GufeTokenizable contains another, the outer object can store the inner
object's key. This allows us to reduce disk usage by only storing one copy
of each object (deduplication).

Deduplication of GufeTokenizables
---------------------------------

There are two types of deduplication of GufeTokenizables. Objects are
deduplicated on storage to disk because we store by reference to the gufe
key. Additionally, objects are deduplicated in memory because we keep a
registry of all instantiated GufeTokenizables.

.. TODO: explain weird sides of that, like the fact that you get the same
   object back when you create another one
