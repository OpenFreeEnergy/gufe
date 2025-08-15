# adapted from from https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/_annotations.py

from typing import Annotated, Any, Type

from openff.units import Quantity
from pydantic import (
    AfterValidator,
    BeforeValidator,
    GetCoreSchemaHandler,
)
from pydantic_core import core_schema

from ..vendor.openff.interchange._annotations import _BoxQuantity as BoxQuantity
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


def make_custom_quantity(unit_name: str) -> Type:
    """Helper function for generating custom quantity types.

    Parameters
    ----------
    unit_name : str
        unit name to validate against (e.g. 'nanometer')

    Returns
    -------
    Type
        A custom type that inherits from openff.units.Quantity.
    """

    CustomQuantity = Annotated[GufeQuantity, AfterValidator(_unit_validator_factory(unit_name))]
    return CustomQuantity


# brute-force these custom types so that mypy recognizes them
NanometerQuantity = Annotated[
    GufeQuantity,
    AfterValidator(_unit_validator_factory("nanometer")),
]

AtmQuantity = Annotated[
    GufeQuantity,
    AfterValidator(_unit_validator_factory("atm")),
]

KelvinQuantity = Annotated[
    GufeQuantity,
    AfterValidator(_unit_validator_factory("kelvin")),
]

# types used elsewhere in the ecosystem
# TODO: add tests here or let that happen in openfe?
NanosecondQuantity = Annotated[
    GufeQuantity,
    AfterValidator(_unit_validator_factory("nanosecond")),
]

PicosecondQuantity = Annotated[
    GufeQuantity,
    AfterValidator(_unit_validator_factory("picosecond")),
]

AngstromQuantity = Annotated[
    GufeQuantity,
    AfterValidator(_unit_validator_factory("angstrom")),
]

KCalPerMolQuantity = Annotated[
    GufeQuantity,
    AfterValidator(_unit_validator_factory("kilocalorie_per_mole")),
]

ArrayQuantity = Annotated[
    GufeQuantity,
    BeforeValidator(_unwrap_list_of_openmm_quantities),
]

NanometerArrayQuantity = Annotated[
    ArrayQuantity,
    AfterValidator(_unit_validator_factory("nanometer")),
]
