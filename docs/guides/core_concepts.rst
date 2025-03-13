Core Concepts
=============

GufeTokenizables
----------------
.. use serialization.rst as starting point

- almost all components (and protocols?) are :class:`.GufeTokenizable`, 
- GufeTokenizables, by definition, are:
  * `immutable <https://docs.python.org/3/glossary.html#term-immutable>`_
  * serializable
  * `hashable <https://docs.python.org/3/glossary.html#term-hashable>`_

Keys
----
.. use serialization.rst as starting point

Flyweight
~~~~~~~~~

Representation
--------------

Any GufeTokenizable can represented in the following ways:

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


Serialization
-------------
.. use serialization.rst as starting point




