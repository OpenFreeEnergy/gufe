# adapted from from https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/_annotations.py

# from enum import StrEnum
from typing import Annotated

from openff.units import Quantity
from pydantic import (
    AfterValidator,
    BeforeValidator,
    WrapSerializer,
    WrapValidator,
)

from ..vendor.openff.interchange._annotations import _BoxQuantity as BoxQuantity
from ..vendor.openff.interchange._annotations import (
    _duck_to_nanometer,
    _unit_validator_factory,
    _unwrap_list_of_openmm_quantities,
    quantity_json_serializer,
    quantity_validator,
)

PydanticQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    WrapSerializer(quantity_json_serializer),
]
AngstromQuantity = Annotated[
    PydanticQuantity,
    AfterValidator(_unit_validator_factory("angstrom")),
]

NanometerQuantity = Annotated[
    PydanticQuantity,
    AfterValidator(_unit_validator_factory("nanometer")),
]

FemtosecondQuantity = Annotated[
    PydanticQuantity,
    AfterValidator(_unit_validator_factory("femtosecond")),
]

NanosecondQuantity = Annotated[
    PydanticQuantity,
    AfterValidator(_unit_validator_factory("nanosecond")),
]

PicosecondQuantity = Annotated[
    PydanticQuantity,
    AfterValidator(_unit_validator_factory("picosecond")),
]

InversePicosecondQuantity = Annotated[
    PydanticQuantity,
    AfterValidator(_unit_validator_factory("1/picosecond")),
]

AtmQuantity = Annotated[
    PydanticQuantity,
    AfterValidator(_unit_validator_factory("atm")),
]

KelvinQuantity = Annotated[
    PydanticQuantity,
    AfterValidator(_unit_validator_factory("kelvin")),
]

KCalPerMolQuantity = Annotated[
    PydanticQuantity,
    AfterValidator(_unit_validator_factory("kilocalorie_per_mole")),
]

TimestepQuantity = Annotated[
    PydanticQuantity,
    AfterValidator(_unit_validator_factory("timestep")),
]

SpringConstantLinearQuantity = Annotated[
    PydanticQuantity,
    AfterValidator(_unit_validator_factory("kilojoule_per_mole / nm ** 2")),
]

SpringConstantAngularQuantity = Annotated[
    PydanticQuantity,
    AfterValidator(_unit_validator_factory("kilojoule_per_mole / radians ** 2")),
]

RadiansQuantity = Annotated[
    PydanticQuantity,
    AfterValidator(_unit_validator_factory("radians")),
]

ArrayQuantity = Annotated[
    PydanticQuantity,
    BeforeValidator(_unwrap_list_of_openmm_quantities),
]

NanometerArrayQuantity = Annotated[
    PydanticQuantity,
    AfterValidator(_unit_validator_factory("nanometer")),
    BeforeValidator(_duck_to_nanometer),
    BeforeValidator(_unwrap_list_of_openmm_quantities),
]

# class CaseInsensitiveStrEnum(StrEnum):
#     # SEE: https://docs.python.org/3/library/enum.html#enum.Enum._missing_
#     @classmethod
#     def _missing_(cls, value):
#         value = value.lower()
#         for member in cls:
#             if member.value == value:
#                 return member
#         return None
