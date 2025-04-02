
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

For our benzene object, that means that its key is directly determined from all items in it's ``to_dict()``
representation, except for ``:version:``, since that is a default parameter:

.. _benzene_to_dict:

.. code:: python

    >>> benzene.defaults()
    {'name': '', ':version:': 1}

    >>> benzene.to_dict()
    {'atoms': [(6, 0, 0, True, 0, 0, {}, 3),
    (6, 0, 0, True, 0, 0, {}, 3),
    (6, 0, 0, True, 0, 0, {}, 3),
    (6, 0, 0, True, 0, 0, {}, 3),
    (6, 0, 0, True, 0, 0, {}, 3),
    (6, 0, 0, True, 0, 0, {}, 3),
    (1, 0, 0, False, 0, 0, {}, 1),
    (1, 0, 0, False, 0, 0, {}, 1),
    (1, 0, 0, False, 0, 0, {}, 1),
    (1, 0, 0, False, 0, 0, {}, 1),
    (1, 0, 0, False, 0, 0, {}, 1),
    (1, 0, 0, False, 0, 0, {}, 1)],
    'bonds': [(0, 1, 12, 0, {}),
    (0, 5, 12, 0, {}),
    (0, 6, 1, 0, {}),
    (1, 2, 12, 0, {}),
    (1, 7, 1, 0, {}),
    (2, 3, 12, 0, {}),
    (2, 8, 1, 0, {}),
    (3, 4, 12, 0, {}),
    (3, 9, 1, 0, {}),
    (4, 5, 12, 0, {}),
    (4, 10, 1, 0, {}),
    (5, 11, 1, 0, {})],
    'conformer': ("\x93NUMPY\x01\x00v\x00{'descr': '<f8', 'fortran_order': False, 'shape': (12, 3), }                                                         \nî|?5^ú9@\x02+\x87\x16ÙN\x15@\x04V\x0e-²\x1d\x13@\x85ëQ¸\x1ee:@²\x9dï§ÆK\x14@Ë¡E¶óý\x0b@×£p=\nW;@q=\n×£p\x17@\x9eï§ÆK7\x07@\x83ÀÊ¡EÖ;@Év¾\x9f\x1a¯\x1b@Zd;ßO\x8d\x0c@ìQ¸\x1e\x85k;@b\x10X9´È\x1c@\x06\x81\x95C\x8bl\x13@sh\x91í|\x7f:@j¼t\x93\x18\x84\x19@ÇK7\x89Aà\x15@í\x9e<,Ô:9@<NÑ\x91\\¾\x12@\x97ÿ\x90~ûú\x14@\x0f\x9c3¢´÷9@\x8d(í\r¾ð\x10@ð\x16HPü\x98\x07@ªñÒMb°;@¼\x05\x12\x14?\x86\x16@Ãdª`TRþ?¦\x9bÄ °\x92<@Ý$\x06\x81\x95C\x1e@Kê\x044\x11¶\x08@RI\x9d\x80&Ò;@\x02\x9a\x08\x1b\x9e\x1e @zÇ):\x92\x8b\x15@9EGrù/:@}?5^ºI\x1a@]mÅþ²û\x19@",
    {}),
    'molprops': {'ofe-name': 'benzene'},
    '__qualname__': 'SmallMoleculeComponent',
    '__module__': 'gufe.components.smallmoleculecomponent',
    ':version:': 1}




This gives the gufe key the following important properties:

* A key is based on a **cryptographic hash**, so it is extremely unlikely
  that two objects that are functionally different will have the same key.
* Key creation is **deterministic**, so that it is preserved across different creation times,
  including across different hardware, across different Python sessions,
  and even within the same Python session.
* A key is preserved across minor versions of the code, since it is dependent on non-default attributes and we follow `SemVer <https://semver.org>`_.

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

.. code:: python

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

1. Serialized Representations of ``GufeTokenizables``
-----------------------------------------------------

- each subclass's implementation of `to_dict()` defines what information gufe will serialize. all other


Representations
^^^^^^^^^^^^^^^

Any GufeTokenizable can be deserialized and

a) dictionary
~~~~~~~~~~~~~

The ``to_dict()`` method is the most explicit way to represent a GufeTokenizable.
This method recursively unpacks any inner GufeTokenizables that an
outer GufeTokenizable contains to their full dict representation.
Although this method is best way to see all information stored in a GufeTokenizable,
it is also the least space-efficient.

For example, we can easily comprehend the ``to_dict()`` representation of benzene :ref:`as shown above <benzene_to_dict>`, but for
a larger and deeply nested object, such as an ``AlchemicalNetwork``, the ``to_dict()`` representation is neither easily readable by humans or computationally memory-efficient.


.. TODO: show this method
.. TODO: diagram

b) shallow dictionary
~~~~~~~~~~~~~~~~~~~~~

The ``to_shallow_dict()`` method is similar to ``to_dict()`` in that it unpacks a tokenizable into a ``dict`` format,
but a shallow dict is *not recursive* and only unpacks the top level of the GufeTokenizable. Anything nested deeper is represented by
the inner objects' GufeTokenizable.

.. code:: python

    >>> alchemical_network.to_shallow_dict()
    {
    'nodes': [
        ChemicalSystem(name=benzene-solvent, components={'ligand': SmallMoleculeComponent(name=benzene), 'solvent': SolventComponent(name=O, K+, Cl-)}),
        ChemicalSystem(name=toluene-solvent, components={'ligand': SmallMoleculeComponent(name=toluene), 'solvent': SolventComponent(name=O, K+, Cl-)}),
        ChemicalSystem(name=styrene-solvent, components={'ligand': SmallMoleculeComponent(name=styrene), 'solvent': SolventComponent(name=O, K+, Cl-)}),
        ChemicalSystem(name=phenol-solvent, components={'ligand': SmallMoleculeComponent(name=phenol), 'solvent': SolventComponent(name=O, K+, Cl-)})
        ],
    'edges': [
        Transformation(stateA=ChemicalSystem(name=benzene-solvent, components={'ligand': SmallMoleculeComponent(name=benzene), 'solvent': SolventComponent(name=O, K+, Cl-)}), stateB=ChemicalSystem(name=toluene-solvent, components={'ligand': SmallMoleculeComponent(name=toluene), 'solvent': SolventComponent(name=O, K+, Cl-)}), protocol=<Protocol-d01baed9cf2500c393bd6ddb35ee38aa>, name=None),
        Transformation(stateA=ChemicalSystem(name=benzene-solvent, components={'ligand': SmallMoleculeComponent(name=benzene), 'solvent': SolventComponent(name=O, K+, Cl-)}), stateB=ChemicalSystem(name=styrene-solvent, components={'ligand': SmallMoleculeComponent(name=styrene), 'solvent': SolventComponent(name=O, K+, Cl-)}), protocol=<Protocol-d01baed9cf2500c393bd6ddb35ee38aa>, name=None),
        Transformation(stateA=ChemicalSystem(name=benzene-solvent, components={'ligand': SmallMoleculeComponent(name=benzene), 'solvent': SolventComponent(name=O, K+, Cl-)}), stateB=ChemicalSystem(name=phenol-solvent, components={'ligand': SmallMoleculeComponent(name=phenol), 'solvent': SolventComponent(name=O, K+, Cl-)}), protocol=<Protocol-d01baed9cf2500c393bd6ddb35ee38aa>, name=None)
        ],
    'name': None,
    '__qualname__': 'AlchemicalNetwork',
    '__module__': 'gufe.network',
    ':version:': 1
    }

.. TODO: diagram


This method is most useful for iterating through the hierarchy of a GufeTokenizable one layer at a time.


c) keyed dictionary
~~~~~~~~~~~~~~~~~~~

The ``to_keyed_dict()`` method is similar to ``to_shallow_dict`` in that it only unpacks the first layer of a GufeTokenizable.
However, a keyed dict represents the next layer as its gufe key, e.g. ``{':gufe-key:': 'ChemicalSystem-96f686efdc070e01b74888cbb830f720'},``

A keyed dict is the most compact representation of a GufeTokenizable and can be useful for understanding its contents,
but it does not have the complete representation for reconstruction or sending information (for this, see the next section, :ref:`keyed chain <keyed_chain>`)

.. code:: python

    >>> alchemical_network.to_keyed_dict()
    {
    'nodes': [
        {':gufe-key:': 'ChemicalSystem-3c648332ff8dccc03a1e1a3d44bc9755'},
        {':gufe-key:': 'ChemicalSystem-655f4d0008a537fe811b11a2dc4a029e'},
        {':gufe-key:': 'ChemicalSystem-6a13159b10c95cb05f542de64ec91fe7'},
        {':gufe-key:': 'ChemicalSystem-ba83a53f18700b3738680da051ff35f3'}
        ],
    'edges': [
        {':gufe-key:': 'Transformation-4d0f802817071c8d14b37efd35187318'},
        {':gufe-key:': 'Transformation-7e7433a86239a41490da52222bf6f78f'},
        {':gufe-key:': 'Transformation-e8d1ccf53116e210d1ccbc3870007271'}
        ],
    'name': None,
    '__qualname__': 'AlchemicalNetwork',
    '__module__': 'gufe.network',
    ':version:': 1
    }


.. TODO: diagram

.. _keyed_chain:

d) keyed chain
~~~~~~~~~~~~~~

The ``to_keyed_chain()`` method is a powerful representation of a GufeTokenizable that enables efficient reconstruction of an object without duplication.
It uses ``to_keyed_dict()`` to unpack a GufeTokenizable from the bottom (innermost) layer up into a flat list of tuples, in the form ``[(gufe_key, keyed_dict)]``. The length of this list is equal to the number of unique GufeTokenizables required to represent the object. This bottom-up deduplication strategy effectively constructs a DAG
(`directed acyclic graph <https://en.wikipedia.org/wiki/Directed_acyclic_graph>`_) where re-used GufeTokenizables are deduplicated.


As an exercise with our example alchemical network, we can look at the first element of each tuple in the keyed dict to see which gufe keys are contained in the alchemical network:

.. code:: python

    >>> [x[0] for x in alchemical_network.to_keyed_chain()()]
    [
    'SolventComponent-e0e47f56b43717156128ad4ae2d49897',
    'SmallMoleculeComponent-3b51f5f92521c712049da092ab061930',
    'SmallMoleculeComponent-ec3c7a92771f8872dab1a9fc4911c795',
    'SmallMoleculeComponent-8225dfb11f2e8157a3fcdcd673d3d40e',
    'Protocol-d01baed9cf2500c393bd6ddb35ee38aa',
    'ChemicalSystem-ba83a53f18700b3738680da051ff35f3',
    'ChemicalSystem-3c648332ff8dccc03a1e1a3d44bc9755',
    'ChemicalSystem-655f4d0008a537fe811b11a2dc4a029e',
    'Transformation-e8d1ccf53116e210d1ccbc3870007271',
    'Transformation-4d0f802817071c8d14b37efd35187318',
    'AlchemicalNetwork-f8bfd63bc848672aa52b081b4d68fadf'
    ]

For keyed chains, the order of the elements in this list matters! When deserializing the 

.. mermaid::

    flowchart TD
        ChemicalSystem-ba83 --> SolventComponent-e0e4
        ChemicalSystem-ba83 --> SmallMoleculeComponent-3
        ChemicalSystem-3c64 --> SolventComponent-e0e4
        ChemicalSystem-3c64 --> SmallMoleculeComponent-e
        ChemicalSystem-655f --> SolventComponent-e0e4
        ChemicalSystem-655f --> SmallMoleculeComponent-8
        Transformation-e8d1 --> ChemicalSystem-3c64
        Transformation-e8d1 --> ChemicalSystem-ba83
        Transformation-e8d1 --> Protocol-d01b

        Transformation-4d0f --> ChemicalSystem-3c64
        Transformation-4d0f --> ChemicalSystem-655f
        Transformation-4d0f --> Protocol-d01b

        AlchemicalNetwork-f8bf --> ChemicalSystem-3c64
        AlchemicalNetwork-f8bf --> ChemicalSystem-655f
        AlchemicalNetwork-f8bf --> ChemicalSystem-ba83
        AlchemicalNetwork-f8bf --> Transformation-4d0f
        AlchemicalNetwork-f8bf --> Transformation-e8d1
        

.. "SolventComponent-e0e4"
.. "SmallMoleculeComponent-3b51"
.. "SmallMoleculeComponent-ec3c"
.. "SmallMoleculeComponent-8225"
.. "Protocol-d01b"
.. "ChemicalSystem-ba83"
.. "ChemicalSystem-3c64"
.. "ChemicalSystem-655f"
.. "Transformation-e8d1"
.. "Transformation-4d0f"
.. "AlchemicalNetwork-f8bf"
.. TODO: maybe show output, maybe abbreviated?
.. TODO: diagram (especially this one!!)



Serialization Methods
^^^^^^^^^^^^^^^^^^^^^

.. TODO explain custom serialization schemes?

- helper methods `to_json` and `to_msgpack` are available, but are simply wrappers around


.. NOTE::
  See :doc:`../how-tos/serialization` for details on how to implement serialization of your own GufeTokenizables.