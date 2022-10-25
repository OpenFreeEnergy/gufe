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
``_defaults``.

Supporting custom types in serialization
----------------------------------------

Occasionally, you may need to add support for a custom type that is not
already known by GUFE. For example, there may be objects that are included
in your GUFE object that come from another library, and are not GUFE
tokenizables.

When this happens, you will need to create and register a custom JSON
codec. 

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

When your object has recursive references
-----------------------------------------

In some cases, your object may have recursive references to other objects.
For example, you may have objects ``parent`` and ``child``, where
``parent.children`` includes ``child`` and ``child.parent`` is ``parent``.
This means that the object dependency graph is not a directed acyclic graph.
The best solution here is to avoid this design pattern whenever possible.
However, if it is unavoidable, 

Understanding the theory: The GUFE key
--------------------------------------

One of the important concepts is that every GUFE object has a unique
identifier, which we call its ``key``. The ``key`` is a string, typically
in the format ``{CLASS_NAME}-{HEXADECIMAL_LABEL}``. The hexadecimal label is
generated based on the contents of the class -- in particular, it is based
on the non-default contents of the ``to_dict`` method. 


Things storage systems may do
-----------------------------

Details on the implementation of storage systems is largely left to
individual projects that build on GUFE. However, the existence of the GUFE
key enables a number of powerful capabilities in storage.
