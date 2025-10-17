# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
# adapted from from https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/_annotations.py
"""
Custom types that inherit from openff.units.Quantity and are pydantic-compatible.
"""

from typing import Annotated, Any, Dict, TypeAlias

import numpy
from openff.units import Quantity
from pydantic import (
    AfterValidator,
    Field,
    PlainSerializer,
    PlainValidator,
    WithJsonSchema,
)

from ..vendor.openff.interchange._annotations import _is_box_shape, _unit_validator_factory


def _plain_quantity_validator(
    value: Any,
) -> Quantity:
    """Take Quantity-like objects and convert them to Quantity objects."""
    # logic adapted from https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/_annotations.py

    def dict_to_quantity(quantity_dict: dict) -> Quantity:
        try:
            # this is coupled to how a Quantity looks in JSON
            return Quantity(quantity_dict["val"], quantity_dict["unit"])
        except KeyError:
            raise ValueError("Quantity must be a dict with keys 'val' and 'unit'.")

    if isinstance(value, Quantity):
        return value
    elif isinstance(value, str):
        return Quantity(value)
    elif isinstance(value, dict):
        return dict_to_quantity(value)
    else:
        raise ValueError(f"Invalid type {type(value)} for Quantity")


def _plain_quantity_serializer(quantity: Quantity) -> Dict[str, Any]:
    return {
        "val": quantity.m,
        "unit": str(quantity.units),
    }


def _is_array(quantity) -> Quantity:
    if isinstance(quantity.m, numpy.ndarray):
        return quantity
    else:
        raise ValueError(f"Quantity {quantity} is not an array.")


# TODO: enforce that this is *not* a list?
GufeQuantity = Annotated[
    Quantity,
    Field(validate_default=True),  # fail fast up front if the default isn't valid.
    PlainValidator(_plain_quantity_validator),
    WithJsonSchema({"type": "number"}),  # this keeps backward compatibility for the JSON schema
    PlainSerializer(_plain_quantity_serializer),
]
"""Pydantic type inherits from ``openff.units.Quantity`` but serializes in a gufe-compatible way."""


GufeArrayQuantity: TypeAlias = Annotated[GufeQuantity, AfterValidator(_is_array)]
"""Convert to a pint.Quantity array, if possible."""


def specify_quantity_units(unit_name: str) -> AfterValidator:
    """Helper function for generating custom quantity types.

    Parameters
    ----------
    unit_name : str
        unit name to validate against (e.g. 'nanometer')

    Returns
    -------
    AfterValidator
        An AfterValidator for defining a custom Quantity type.

    """

    return AfterValidator(_unit_validator_factory(unit_name))


NanometerQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("nanometer"),
]
"""Pydantic type that requires a ``pint.Quantity`` compatible with nanometers. Input will be converted to nanometers upon model validation."""

BarQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("bar"),
]
"""Pydantic type that requires a ``pint.Quantity`` compatible with bar. Input will be converted to bar upon model validation."""

KelvinQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("kelvin"),
]
"""Pydantic type that requires a ``pint.Quantity`` compatible with kelvin. Input will be converted to kelvin upon model validation.
Note: to define input in Celsius, you must use ``Quantity`` explicitly, e.g. ``openff.units.Quantity(25, "celsius")`` instead of ``25 * unit.celsius``.
See  https://pint.readthedocs.io/en/stable/user/nonmult.html for more information.
"""

# types used elsewhere in the ecosystem
NanosecondQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("nanosecond"),
]
"""Pydantic type that requires a ``pint.Quantity`` compatible with nanoseconds. Input will be converted to nanoseconds upon model validation."""

PicosecondQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("picosecond"),
]
"""Pydantic type that requires a ``pint.Quantity`` compatible with picoseconds. Input will be converted to picoseconds upon model validation."""

AngstromQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("angstrom"),
]
"""Pydantic type that requires a ``pint.Quantity`` compatible with angstroms. Input will be converted to angstroms upon model validation."""

KCalPerMolQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("kilocalorie_per_mole"),
]
"""Pydantic type that requires a ``pint.Quantity`` compatible with kilocalorie_per_mole. Input will be converted to kilocalorie_per_mole upon model validation."""

VoltsQuantity: TypeAlias = Annotated[GufeQuantity, specify_quantity_units("volts")]
"""Pydantic type that requires a ``pint.Quantity`` compatible with volts. Input will be converted to volts upon model validation."""

NanometerArrayQuantity: TypeAlias = Annotated[
    GufeArrayQuantity,
    specify_quantity_units("nanometer"),
]
"""Pydantic type that requires an array of ``pint.Quantity`` compatible with nanometer. Input will be converted to nanometers upon model validation."""

BoxQuantity = Annotated[
    NanometerQuantity,
    AfterValidator(_is_box_shape),
]
"""
Pydantic type that requires a 3x3 or 3x1 array of ``pint.Quantity`` values compatible with nanometers.
Input will be converted to nanometer upon model validation.
If input is a 3x1 array, it will be converted to 3x3 diagonal array with zeroes on the off-diagonal.
"""
