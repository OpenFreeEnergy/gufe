# adapted from from https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/_annotations.py
"""
Custom types that inherit from openff.units.Quantity and are pydantic-compatible.
"""

from typing import Annotated, Any, TypeAlias

from openff.units import Quantity
from pydantic import (
    AfterValidator,
    BeforeValidator,
    GetCoreSchemaHandler,
)
from pydantic_core import core_schema

from ..vendor.openff.interchange._annotations import (
    _unit_validator_factory,
    _unwrap_list_of_openmm_quantities,
    quantity_json_serializer,
    quantity_validator,
)


class _QuantityPydanticAnnotation:
    @classmethod
    def __get_pydantic_core_schema__(
        cls,
        source: Any,
        handler: GetCoreSchemaHandler,
    ) -> core_schema.CoreSchema:
        """
        This Annotation lets us define a GufeQuantity that is identical to
        an openff-units Quantity, except it's also pydantic-compatible.
        """
        json_schema = core_schema.with_info_wrap_validator_function(
            function=quantity_validator, schema=core_schema.float_schema()
        )
        python_schema = core_schema.with_info_wrap_validator_function(
            function=quantity_validator,
            schema=core_schema.is_instance_schema(Quantity),
        )
        serialize_schema = core_schema.wrap_serializer_function_ser_schema(quantity_json_serializer)
        return core_schema.json_or_python_schema(
            json_schema=json_schema,
            python_schema=python_schema,
            serialization=serialize_schema,
        )


GufeQuantity = Annotated[Quantity, _QuantityPydanticAnnotation]


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
"""Convert a pint.Quantity or to nanometers, if possible."""

AtmQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("atm"),
]
"""Convert a pint.Quantity or to atm, if possible."""

KelvinQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("kelvin"),
]
"""Convert a pint.Quantity or to kelvin, if possible."""

# types used elsewhere in the ecosystem
NanosecondQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("nanosecond"),
]
"""Convert a pint.Quantity or to nanoseconds, if possible."""


PicosecondQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("picosecond"),
]
"""Convert a pint.Quantity or to picoseconds, if possible."""

AngstromQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("angstrom"),
]
"""Convert a pint.Quantity or to angstroms, if possible."""

KCalPerMolQuantity: TypeAlias = Annotated[
    GufeQuantity,
    specify_quantity_units("kilocalorie_per_mole"),
]
"""Convert a pint.Quantity or to kcal/mol, if possible."""

GufeArrayQuantity: TypeAlias = Annotated[
    GufeQuantity,
    BeforeValidator(_unwrap_list_of_openmm_quantities),
]

NanometerArrayQuantity: TypeAlias = Annotated[
    GufeArrayQuantity,
    specify_quantity_units("nanometer"),
]
