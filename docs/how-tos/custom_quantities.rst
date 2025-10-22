.. _howto-quantity:

How to define a custom ``GufeQuantity``
=======================================

:class:`gufe.settings.typing.GufeQuantity` is compatible with both pydantic serialization and openff.

See :class:`gufe.settings.typing` for quantities included in gufe.

To make your own:


.. code:: python

    from gufe.settings.typing import GufeQuantity, specify_quantity_units
    RadiansQuantity:TypeAlias = Annotated[GufeQuantity, specify_quantity_units("radians")]


or for an array:

.. code::  python

    from gufe.settings.typing import GufeArrayQuantity, specify_quantity_units

    RadiansArrayQuantity: TypeAlias = Annotated[
        GufeArrayQuantity,
        specify_quantity_units("radians"),
    ]
