# adapted from from https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/_annotations.py

from enum import StrEnum
from typing import Annotated

from openff.toolkit import Quantity
from pydantic import (
    AfterValidator,
    WrapSerializer,
    WrapValidator,
)

from ..vendor.openff.interchange._annotations import (
    _unit_validator_factory,
    quantity_json_serializer,
    quantity_validator,
)

NanometerQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    AfterValidator(_unit_validator_factory("nanometer")),
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

# TODO: openfe will need this.
# _AngstromQuantity = Annotated[
#     Quantity,
#     WrapValidator(quantity_validator),
#     AfterValidator(_is_angstrom),
#     WrapSerializer(quantity_json_serializer),
# ]


# class CaseInsensitiveStrEnum(StrEnum):
#     # SEE: https://docs.python.org/3/library/enum.html#enum.Enum._missing_
#     @classmethod
#     def _missing_(cls, value):
#         value = value.lower()
#         for member in cls:
#             if member.value == value:
#                 return member
#         return None
