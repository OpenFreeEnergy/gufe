"""Custom models for dealing with unit-bearing quantities in a Pydantic-compatible manner."""

import json
from typing import TYPE_CHECKING, Any

import numpy
from openff.units import Quantity, Unit, unit
from openff.utilities import has_package, requires_package

from .exceptions import (
    MissingUnitError,
    UnitValidationError,
    UnsupportedExportError,
)

if TYPE_CHECKING:
    import openmm.unit


class _FloatQuantityMeta(type):
    def __getitem__(self, t):
        return type("FloatQuantity", (FloatQuantity,), {"__unit__": t})


if TYPE_CHECKING:
    FloatQuantity = unit.Quantity
else:

    class FloatQuantity(float, metaclass=_FloatQuantityMeta):
        """A model for unit-bearing floats."""

        @classmethod
        def __get_validators__(cls):
            yield cls.validate_type

        @classmethod
        def validate_type(cls, val):
            """Process a value tagged with units into one tagged with "OpenFF" style units."""
            unit_ = getattr(cls, "__unit__", Any)
            if unit_ is Any:
                if isinstance(val, (float, int)):
                    # TODO: Can this exception be raised with knowledge of the field it's in?
                    raise MissingUnitError(f"Value {val} needs to be tagged with a unit")
                elif isinstance(val, Quantity):
                    return Quantity(val)
                elif _is_openmm_quantity(val):
                    return _from_omm_quantity(val)
                else:
                    raise UnitValidationError(f"Could not validate data of type {type(val)}")
            else:
                unit_ = Unit(unit_)
                if isinstance(val, Quantity):
                    # some custom behavior could go here
                    assert unit_.dimensionality == val.dimensionality
                    # return through converting to some intended default units (taken from the class)
                    val._magnitude = float(val.m)
                    return val.to(unit_)

                if _is_openmm_quantity(val):
                    return _from_omm_quantity(val).to(unit_)
                if isinstance(val, int) and not isinstance(val, bool):
                    # coerce ints into floats for a FloatQuantity
                    return float(val) * unit_
                if isinstance(val, float):
                    return val * unit_
                if isinstance(val, str):
                    # could do custom deserialization here?
                    val = Quantity(val).to(unit_)
                    val._magnitude = float(val._magnitude)
                    return val
                if "unyt" in str(val.__class__):
                    if val.value.shape == ():
                        # this is a scalar force into an array by unyt's design
                        if "float" in str(val.value.dtype):
                            return float(val.value) * unit_
                        elif "int" in str(val.value.dtype):
                            return int(val.value) * unit_

                raise UnitValidationError(f"Could not validate data of type {type(val)}")


def _is_openmm_quantity(obj: object) -> bool:
    if has_package("openmm"):
        import openmm.unit

        return isinstance(obj, openmm.unit.Quantity)

    else:
        return "openmm.unit.quantity.Quantity" in str(type(object))


@requires_package("openmm.unit")
def _from_omm_quantity(val: "openmm.unit.Quantity") -> Quantity:
    """
    Convert float or array quantities tagged with SimTK/OpenMM units to a Pint-compatible quantity.
    """
    unit_: openmm.unit.Unit = val.unit
    val_ = val.value_in_unit(unit_)
    if type(val_) in {float, int}:
        unit_ = val.unit
        return float(val_) * Unit(str(unit_))
    # Here is where the toolkit's ValidatedList could go, if present in the environment
    elif (type(val_) in {tuple, list, numpy.ndarray}) or (type(val_).__module__ == "openmm.vec3"):
        array = numpy.asarray(val_)
        return array * Unit(str(unit_))
    elif isinstance(val_, (float, int)) and type(val_).__module__ == "numpy":
        return val_ * Unit(str(unit_))
    else:
        raise UnitValidationError(
            "Found a openmm.unit.Unit wrapped around something other than a float-like "
            f"or numpy.ndarray-like. Found a unit wrapped around type {type(val_)}."
        )


class QuantityEncoder(json.JSONEncoder):
    """
    JSON encoder for unit-wrapped floats and NumPy arrays.

    This is intended to operate on FloatQuantity and ArrayQuantity objects.
    """

    def default(self, obj):
        if isinstance(obj, Quantity):
            if isinstance(obj.magnitude, (float, int)):
                data = obj.magnitude
            elif isinstance(obj.magnitude, numpy.ndarray):
                data = obj.magnitude.tolist()
            else:
                # This shouldn't ever be hit if our object models
                # behave in ways we expect?
                raise UnsupportedExportError(f"trying to serialize unsupported type {type(obj.magnitude)}")
            return {
                "val": data,
                "unit": str(obj.units),
            }


def custom_quantity_encoder(v):
    """Wrap json.dump to use QuantityEncoder."""
    return json.dumps(v, cls=QuantityEncoder)


def json_loader(data: str) -> dict:
    """Load JSON containing custom unit-tagged quantities."""
    # TODO: recursively call this function for nested models
    out: dict = json.loads(data)
    for key, val in out.items():
        try:
            # Directly look for an encoded FloatQuantity/ArrayQuantity,
            # which is itself a dict
            v = json.loads(val)
        except (json.JSONDecodeError, TypeError):
            # Handles some cases of the val being a primitive type
            continue
        # TODO: More gracefully parse non-FloatQuantity/ArrayQuantity dicts
        unit_ = Unit(v["unit"])
        val = v["val"]
        out[key] = unit_ * val
    return out


class _ArrayQuantityMeta(type):
    def __getitem__(self, t):
        return type("ArrayQuantity", (ArrayQuantity,), {"__unit__": t})


if TYPE_CHECKING:
    ArrayQuantity = unit.Quantity
else:

    class ArrayQuantity(float, metaclass=_ArrayQuantityMeta):
        """A model for unit-bearing arrays."""

        @classmethod
        def __get_validators__(cls):
            yield cls.validate_type

        @classmethod
        def validate_type(cls, val):
            """Process an array tagged with units into one tagged with "OpenFF" style units."""
            unit_ = getattr(cls, "__unit__", Any)
            if unit_ is Any:
                if isinstance(val, (list, numpy.ndarray)):
                    # Work around a special case in which val might be list[openmm.unit.Quantity]
                    if isinstance(val, list) and {type(element).__module__ for element in val} == {
                        "openmm.unit.quantity"
                    }:
                        unit_ = _from_omm_quantity(val[-1]).units
                        return Quantity(
                            [_from_omm_quantity(element).m for element in val],
                            units=unit_,
                        )

                    # TODO: Can this exception be raised with knowledge of the field it's in?
                    raise MissingUnitError(f"Value {val} needs to be tagged with a unit")

                elif isinstance(val, Quantity):
                    # TODO: This might be a redundant cast causing wasted CPU time.
                    #       But maybe it handles pint vs openff.units.unit?
                    return Quantity(val)
                elif _is_openmm_quantity(val):
                    return _from_omm_quantity(val)
                else:
                    raise UnitValidationError(f"Could not validate data of type {type(val)}")
            else:
                unit_ = Unit(unit_)
                if isinstance(val, Quantity):
                    assert unit_.dimensionality == val.dimensionality
                    return val.to(unit_)
                if _is_openmm_quantity(val):
                    return _from_omm_quantity(val).to(unit_)
                if isinstance(val, (numpy.ndarray, list)):
                    if "unyt" in str(val.__class__):
                        val = val.to_ndarray()
                    try:
                        return val * unit_
                    except RuntimeError as error:
                        # unyt subclasses ndarray but doesn't __mult__ with
                        # pint.Unit objects
                        if val.__class__.__module__.startswith("unyt"):
                            return val.to_ndarray() * unit_
                        else:
                            raise error
                if isinstance(val, bytes):
                    # Define outside loop
                    dt = numpy.dtype(int).newbyteorder("<")
                    return numpy.frombuffer(val, dtype=dt) * unit_
                if isinstance(val, str):
                    # could do custom deserialization here?
                    raise NotImplementedError
                raise UnitValidationError(f"Could not validate data of type {type(val)}")
