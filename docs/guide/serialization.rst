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
  :class:`.JSONCodec` has been registered (see :ref:`Supporting custom types
  in serialization`).
* ``_from_dict`` is a ``classmethod`` that takes the dictionary returned by
  ``_to_dict`` and reconstructs the object.
* ``_defaults`` returns whatever ``_to_dict`` would return for any default
  values; i.e., it should be a dictionary ...


As a very simple example, consider the following complete implementation:

.. code::

    class Foo(GufeTokenizable):
        def __init__(self, bar, baz=None):
            self.bar = bar
            self.baz = baz

        def _to_dict(self):
            return {'bar': self.bar, 'baz': self.baz}

        @classmethod
        def _from_dict(cls, dct):
            return cls(**dct)

        def _defaults(self):
            return {'baz': None}


Supporting custom types in serialization
----------------------------------------

Occasionally, you may need to add support for a custom type that is not
already known by GUFE. For example, there may be objects that are included
in your GUFE object that come from another library, and are not GUFE
tokenizables.

When this happens, you will need to create and register a custom JSON
codec. To do this, at a minimum you need to define a ``to_dict`` and
``from_dict`` method for your type. Additionally, the 

Dumping arbitrary objects to JSON
---------------------------------

Advanced storage systems are space-efficient by only storing each object
once. However, any :class:`.GufeTokenizable` can be dumped to a JSON file
using the custom JSON handlers. Given a :class:`.GufeTokenizable` called
``obj`` and a path-like called ``filename``, you can dump to JSON with this
recipe:

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

Using JSON codecs outside of JSON
---------------------------------

In a custom recursive storage scheme, 

When your object has recursive references
-----------------------------------------

In some cases, your object may have recursive references to other objects.
For example, you may have objects ``parent`` and ``child``, where
``parent.children`` includes ``child`` and ``child.parent`` is ``parent``.
This means that the object dependency graph is not a directed acyclic graph.
The best solution here is to avoid this design pattern whenever possible.
However, if it is unavoidable, the way to do this is to set the
``child.parent`` attribute as part of the ``parent._from_dict`` method.

Understanding the theory: The GUFE key
--------------------------------------

One of the important concepts is that every GUFE object has a unique
identifier, which we call its ``key``. The ``key`` is a string, typically
in the format ``{CLASS_NAME}-{HEXADECIMAL_LABEL}``. The hexadecimal label is
generated based on the contents of the class -- in particular, it is based
on contents of the ``_to_dict`` dictionary, filtered to remove anything that
matches the ``_defaults`` dictionary.

This gives the GUFE key a number of important properties:

* GUFE keys can be used for store-by-reference: When one GufeTokenizable
  contains another, the inner object can store 

Deduplication of GufeTokenizables
---------------------------------
