Core Concepts
=============

GufeTokenizables
----------------

- almost all components (and protocols?) are :class:`.GufeTokenizable`, 
- GufeTokenizables, by definition, are:
  * `immutable <https://docs.python.org/3/glossary.html#term-immutable>`_
  * serializable
  * `hashable <https://docs.python.org/3/glossary.html#term-hashable>`_


Representation
~~~~~~~~~~~~~~

Any GufeTokenizable can represented in the following ways:

1. dict
  - this is the most explicit way to represent a GufeTokenizable
  - unpacks all levels to their dict representation
2. keyed_dict
  - every GufeTokenizable will be represented as its key
3. shallow_dict
  - only one level is 'unpacked', anything deeper is stored by key 
  - (QUESTION: i thought? in tokenization.py, `to_shallow_dict` is just calling to_dict?)
4. keyed_chain
  - creates a DAG for efficient reconstruction without duplication


Serialization
~~~~~~~~~~~~~


Keys
~~~~


