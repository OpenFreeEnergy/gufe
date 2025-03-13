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
    from gufe.tokenization import JSON_HANDLER, GufeTokenizable
    with open(filename, mode='r') as f:
        obj = GufeTokenizable.from_dict(json.load(f, cls=JSON_HANDLER.decoder))

Note that these objects are not space-efficient: that is, if you have
the same object in memory referenced by multiple objects (e.g., an identical
``ProteinComponent`` in more than one ``ChemicalSystem``), then you will
save multiple copies of its JSON representation.

On reloading, tools that use the recommended ``from_dict`` method will undo
this duplication; see :ref:`gufe-memory-deduplication` for details.

As a more space-efficient alternative to ``to_dict``/``from_dict``, consider
using ``to_keyed_chain``/``from_keyed_chain`` instead.
This deals in a representation using the :class:`.KeyedChain` approach, which
avoids duplication of dependent :class:`.GufeTokenizables` in the serialized
JSON representation.

Convenient serialization
~~~~~~~~~~~~------------

We also provide convenience methods to convert any :class:`.GufeTokenizable` to
and from JSON using a space-efficient serialization strategy based on our
:class:`.KeyedChain` representation. This is intended for developers that want
to serialise these objects using the current best practice and are not
concerned with the details of the process. The :func:`to_json
<gufe.tokenization.GufeTokenizable.to_json>` API offers the flexibility to
convert to JSON directly or to write to a filelike object:

.. code::

    # get a json representation in-memory
    json = obj.to_json()

    # save to a file directly
    obj.to_json(file=filename)

Similarly, you can recreate the object using the :func:`from_json <gufe.tokenization.GufeTokenizable.from_json>`
classmethod:

.. code::

    # load the object from a json file produced with `to_json`
    obj = cls.from_json(file=filename)

    # load from a string produced with `to_json`
    obj = cls.from_json(content=json)

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