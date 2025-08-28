# adapted from from https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/_annotations.py
"""
Custom types that inherit from openff.units.Quantity and are pydantic-compatible.
"""

from typing import Annotated, Any, Dict, TypeAlias

import numpy
from openff.units import Quantity
from pydantic import (
    AfterValidator,
    BeforeValidator,
    Field,
    PlainSerializer,
    PlainValidator,
    ValidationInfo,
    WithJsonSchema,
)

from ..vendor.openff.interchange._annotations import (
    _duck_to_nanometer,
    _is_box_shape,
    _unit_validator_factory,
    _unwrap_list_of_openmm_quantities,
)


def _plain_quantity_validator(
    value: Any,
    info: ValidationInfo,
) -> Quantity:
    """Take Quantity-like objects and convert them to Quantity objects."""

    # logic from https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/_annotations.py
    if info.mode == "json":
        assert isinstance(value, dict), "Quantity must be in dict form here."
        # this is coupled to how a Quantity looks in JSON
        return Quantity(value["val"], value["unit"])

        # some more work may be needed to work with arrays, lists, tuples, etc.

    assert info.mode == "python"

    if isinstance(value, Quantity):
        return value
    elif isinstance(value, str):
        return Quantity(value)
    elif isinstance(value, dict):
        return Quantity(value["val"], value["unit"])
    if "openmm" in str(type(value)):
        from openff.units.openmm import from_openmm

        return from_openmm(value)
    else:
        raise ValueError(f"Invalid type {type(value)} for Quantity")


def _plain_quantity_serializer(quantity: Quantity) -> Dict[str, Any]:
    # logic from https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/_annotations.py
    magnitude = quantity.m

    if isinstance(magnitude, numpy.ndarray):
        # This could be something fancier, list a bytestring
        magnitude = magnitude.tolist()

    return {
        "val": magnitude,
        "unit": str(quantity.units),
    }


GufeQuantity = Annotated[
    Quantity,
    Field(validate_default=True),  # fail fast up front if the default isn't valid.
    PlainValidator(_plain_quantity_validator),
    WithJsonSchema({"type": "number"}),  # this keeps backward compatibility for the JSON schema
    PlainSerializer(_plain_quantity_serializer),
]


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

GufeArrayQuantity: TypeAlias = Annotated[
    GufeQuantity,
    BeforeValidator(_unwrap_list_of_openmm_quantities),
]
"""Convert to a list of pint.Quantity objects, if possible."""

NanometerArrayQuantity: TypeAlias = Annotated[
    GufeArrayQuantity,
    specify_quantity_units("nanometer"),
]
"""Convert to a list of pint.Quantity objects to nanometers, if possible."""

BoxQuantity = Annotated[
    NanometerQuantity,
    AfterValidator(_is_box_shape),
    BeforeValidator(_duck_to_nanometer),
    BeforeValidator(_unwrap_list_of_openmm_quantities),
]
