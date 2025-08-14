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
    _unit_validator_factory as unit_validator,
    _unwrap_list_of_openmm_quantities,
    quantity_json_serializer,
    quantity_validator,
)

GufeQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    WrapSerializer(quantity_json_serializer),
]

NanometerQuantity = Annotated[
    GufeQuantity,
    AfterValidator(unit_validator("nanometer")),
]

AtmQuantity = Annotated[
    GufeQuantity,
    AfterValidator(unit_validator("atm")),
]

KelvinQuantity = Annotated[
    GufeQuantity,
    AfterValidator(unit_validator("kelvin")),
]

# NanosecondQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(unit_validator("nanosecond")),
# ]

# AngstromQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(unit_validator("angstrom")),
# ]

# FemtosecondQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(unit_validator("femtosecond")),
# ]

# PicosecondQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(unit_validator("picosecond")),
# ]

# InversePicosecondQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(unit_validator("1/picosecond")),
# ]

# KCalPerMolQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(unit_validator("kilocalorie_per_mole")),
# ]

# TimestepQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(unit_validator("timestep")),
# ]

# SpringConstantLinearQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(unit_validator("kilojoule_per_mole / nm ** 2")),
# ]

# SpringConstantAngularQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(unit_validator("kilojoule_per_mole / radians ** 2")),
# ]

# RadiansQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(unit_validator("radians")),
# ]

ArrayQuantity = Annotated[
    GufeQuantity,
    BeforeValidator(_unwrap_list_of_openmm_quantities),
]

NanometerArrayQuantity = Annotated[
    ArrayQuantity,
    AfterValidator(unit_validator("nanometer")),
]