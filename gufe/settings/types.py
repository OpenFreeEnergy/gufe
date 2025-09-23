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
"""Convert a pint.Quantity to nanometers, if possible."""

AtmQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("atm"),
]
"""Convert a pint.Quantity to atm, if possible."""

KelvinQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("kelvin"),
]
"""Convert a pint.Quantity to kelvin, if possible."""

# types used elsewhere in the ecosystem
NanosecondQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("nanosecond"),
]
"""Convert a pint.Quantity to nanoseconds, if possible."""

PicosecondQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("picosecond"),
]
"""Convert a pint.Quantity to picoseconds, if possible."""

AngstromQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("angstrom"),
]
"""Convert a pint.Quantity to angstroms, if possible."""

KCalPerMolQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("kilocalorie_per_mole"),
]
"""Convert a pint.Quantity to kcal/mol, if possible."""

VoltsQuantity: TypeAlias = Annotated[GufeQuantity, specify_quantity_units("volts")]
"""Convert a pint.Quantity to volts, if possible."""

NanometerArrayQuantity: TypeAlias = Annotated[
    GufeArrayQuantity,
    specify_quantity_units("nanometer"),
]
"""Convert to a list of pint.Quantity objects to nanometers, if possible."""

BoxQuantity = Annotated[
    NanometerQuantity,
    AfterValidator(_is_box_shape),
]
