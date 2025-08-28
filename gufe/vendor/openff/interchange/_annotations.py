# Vendored from https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/_annotations.py
import functools
from collections.abc import Callable
from typing import Annotated, Any

import numpy
from annotated_types import Gt
from openff.units import Quantity  # import from units so we don't have to build toolkit just for docs

PositiveFloat = Annotated[float, Gt(0)]


def _has_compatible_dimensionality(
    quantity: Quantity,
    unit: str,
    convert: bool,
) -> Quantity:
    """Check if a Quantity has the same dimensionality as a given unit and optionally convert."""
    if quantity.is_compatible_with(unit):
        if convert:
            return quantity.to(unit)
        else:
            return quantity
    else:
        raise ValueError(
            f"Dimensionality of {quantity=} is not compatible with {unit=}",
        )


def _unit_validator_factory(unit: str) -> Callable:
    """Return a function, meant to be passed to a validator, that checks for a specific unit."""
    return functools.partial(_has_compatible_dimensionality, unit=unit, convert=True)


def _duck_to_nanometer(value: Any):
    """Cast list or ndarray without units to Quantity[ndarray] of nanometer."""
    if isinstance(value, (list, numpy.ndarray)):
        return Quantity(value, "nanometer")
    else:
        return value


def _is_box_shape(quantity) -> Quantity:
    if quantity.m.shape == (3, 3):
        return quantity
    elif quantity.m.shape == (3,):
        return numpy.eye(3) * quantity
    else:
        raise ValueError(f"Quantity {quantity} is not a box.")


def _unwrap_list_of_openmm_quantities(value: Any):
    """Unwrap a list of OpenMM quantities to a single Quantity."""
    if isinstance(value, list):
        if any(["openmm" in str(type(element)) for element in value]):
            from openff.units.openmm import from_openmm

            if len({element.unit for element in value}) != 1:
                raise ValueError("All units must be the same.")

            return from_openmm(value)

        else:
            return value

    else:
        return value
