.. _howto-quantity:

How to define a custom ``GufeQuantity``
=====================================

``GufeQuantity`` is compatible with both pydantic serialization and openff.

To make your own:


.. code:: python

    from gufe.settings.types import GufeQuantity, specify_quantity_units
    RadiansQuantity:TypeAlias = Annotated[GufeQuantity, specify_quantity_units("radians")]


or for an array:

.. code:: python
    from gufe.settings.types import GufeArrayQuantity, specify_quantity_units

    RadiansArrayQuantity: TypeAlias = Annotated[
    GufeArrayQuantity,
    specify_quantity_units("radians"),
    ]
