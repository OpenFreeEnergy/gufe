.. _concepts-serialization:

Serialization constraints and methods in **gufe**
=================================================

The ``GufeTokenizable`` object is designed to be easily serialized by a variety of methods (see also: :ref:`concepts-tokenizables-serialization-methods`).
In order these objects to be serializable by these methods, however, there exist constraints on the data types that a ``GufeTokenizable`` can be composed of.

In addition to the data types natively supported by `JSON <https://docs.python.org/3/library/json.html#encoders-and-decoders>`_ and `MessagePack <https://github.com/msgpack/msgpack/blob/master/spec.md#serialization-type-to-format-conversion>`_, **gufe** also supports the following types as serializable attributes:

- ``pathlib.Path``
- ``numpy.generic``
- ``numpy.ndarray``
- ``bytes``
- ``datetime.datetime``
- ``openff.units.Unit``
- ``openff.units.Quantity``
- ``uuid.UUID``
- ``gufe.settings.SettingsBaseModel``
