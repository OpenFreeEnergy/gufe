# adapted from from https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/_annotations.py


from typing import Annotated

from openff.toolkit import Quantity

from gufe.vendor.openff.interchange._annotations import (
    quantity_json_serializer,
    quantity_validator,
    _unit_validator_factory,
)

from pydantic import (
    AfterValidator,
    WrapSerializer,
    WrapValidator,
)


(
    _is_dimensionless,
    _is_kj_mol,
    _is_nanometer,
    _is_atm,
    _is_angstrom,
    _is_degree,
    _is_elementary_charge,
) = (
    _unit_validator_factory(unit=_unit)
    for _unit in [
        "dimensionless",
        "kilojoule / mole",
        "nanometer",
        "atm",
        "angstrom",
        "degree",
        "elementary_charge",
    ]
)

_NanometerQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    AfterValidator(_is_nanometer),
    WrapSerializer(quantity_json_serializer),
]

_PressureQuantity = Annotated[
    Quantity,
    WrapValidator(quantity_validator),
    AfterValidator(_is_atm),
    WrapSerializer(quantity_json_serializer),
]

# _AngstromQuantity = Annotated[
#     Quantity,
#     WrapValidator(quantity_validator),
#     AfterValidator(_is_angstrom),
#     WrapSerializer(quantity_json_serializer),
# ]
