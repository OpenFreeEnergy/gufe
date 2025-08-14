# adapted from from https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/_annotations.py

from typing import Annotated, Type

from openff.units import Quantity
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

GufeQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    WrapSerializer(quantity_json_serializer),
]


def custom_quantity(unit_name: str) -> Type:
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

    CustomQuantity = Annotated[
        GufeQuantity, AfterValidator(_unit_validator_factory(unit_name))
    ]
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

AngstromQuantity = Annotated[
    GufeQuantity,
    AfterValidator(_unit_validator_factory("angstrom")),
]

# FemtosecondQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(_unit_validator_factory("femtosecond")),
# ]

# PicosecondQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(_unit_validator_factory("picosecond")),
# ]

# InversePicosecondQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(_unit_validator_factory("1/picosecond")),
# ]

# KCalPerMolQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(_unit_validator_factory("kilocalorie_per_mole")),
# ]

# TimestepQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(_unit_validator_factory("timestep")),
# ]

# SpringConstantLinearQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(_unit_validator_factory("kilojoule_per_mole / nm ** 2")),
# ]

# SpringConstantAngularQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(_unit_validator_factory("kilojoule_per_mole / radians ** 2")),
# ]

# RadiansQuantity = Annotated[
#     GufeQuantity,
#     AfterValidator(_unit_validator_factory("radians")),
# ]

ArrayQuantity = Annotated[
    GufeQuantity,
    BeforeValidator(_unwrap_list_of_openmm_quantities),
]

NanometerArrayQuantity = Annotated[
    ArrayQuantity,
    AfterValidator(_unit_validator_factory("nanometer")),
]