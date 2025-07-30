# adapted from from https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/_annotations.py

# from enum import StrEnum
from typing import Annotated

from openff.toolkit import Quantity
from pydantic import (
    AfterValidator,
    BeforeValidator,
    WrapSerializer,
    WrapValidator,
)

from ..vendor.openff.interchange._annotations import _BoxQuantity as BoxQuantity
from ..vendor.openff.interchange._annotations import (
    _unit_validator_factory,
    _unwrap_list_of_openmm_quantities,
    quantity_json_serializer,
    quantity_validator,
)

AngstromQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    AfterValidator(_unit_validator_factory("angstrom")),
    WrapSerializer(quantity_json_serializer),
]

NanometerQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    AfterValidator(_unit_validator_factory("nanometer")),
    WrapSerializer(quantity_json_serializer),
]

FemtosecondQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    AfterValidator(_unit_validator_factory("femtosecond")),
    WrapSerializer(quantity_json_serializer),
]

NanosecondQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    AfterValidator(_unit_validator_factory("nanosecond")),
    WrapSerializer(quantity_json_serializer),
]

PicosecondQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    AfterValidator(_unit_validator_factory("picosecond")),
    WrapSerializer(quantity_json_serializer),
]

InversePicosecondQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    AfterValidator(_unit_validator_factory("1/picosecond")),
    WrapSerializer(quantity_json_serializer),
]

AtmQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    AfterValidator(_unit_validator_factory("atm")),
    WrapSerializer(quantity_json_serializer),
]

KelvinQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    AfterValidator(_unit_validator_factory("kelvin")),
    WrapSerializer(quantity_json_serializer),
]

KCalPerMolQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    AfterValidator(_unit_validator_factory("kilocalorie_per_mole")),
    WrapSerializer(quantity_json_serializer),
]

TimestepQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    AfterValidator(_unit_validator_factory("timestep")),
    WrapSerializer(quantity_json_serializer),
]

SpringConstantLinearQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    AfterValidator(_unit_validator_factory("kilojoule_per_mole / nm ** 2")),
    WrapSerializer(quantity_json_serializer),
]

SpringConstantAngularQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    AfterValidator(_unit_validator_factory("kilojoule_per_mole / radians ** 2")),
    WrapSerializer(quantity_json_serializer),
]

RadiansQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    AfterValidator(_unit_validator_factory("radians")),
    WrapSerializer(quantity_json_serializer),
]

ArrayQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    BeforeValidator(_unwrap_list_of_openmm_quantities),
    WrapSerializer(quantity_json_serializer),
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
