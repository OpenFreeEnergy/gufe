# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""
Pydantic models used for storing settings.
"""

import abc
import pprint
from typing import Annotated, Any, Literal, TypeAlias

from annotated_types import Ge
from openff.units import unit
from pydantic import BeforeValidator, ConfigDict, Field, InstanceOf, PositiveFloat, PrivateAttr, field_validator

from ..vendor.openff.interchange.pydantic import _BaseModel
from .types import AtmQuantity, GufeQuantity, KelvinQuantity, NanometerQuantity, specify_quantity_units

VoltsQuantity: TypeAlias = Annotated[GufeQuantity, specify_quantity_units("volts")]

class SettingsBaseModel(_BaseModel):
    """Settings and modifications we want for all settings classes."""

    _is_frozen: bool = PrivateAttr(default_factory=lambda: False)
    model_config = ConfigDict(
        extra="forbid",
        arbitrary_types_allowed=True,  # needed to parse custom types
    )

    def _ipython_display_(self):
        pprint.pprint(self.dict())

    def frozen_copy(self):
        """A copy of this Settings object which cannot be modified

        This is intended to be used by Protocols to make their stored Settings
        read-only
        """
        copied = self.model_copy(deep=True)

        def freeze_model(model):
            submodels = (
                mod for field in model.model_fields if isinstance(mod := getattr(model, field), SettingsBaseModel)
            )
            for mod in submodels:
                freeze_model(mod)

            if not model._is_frozen:
                model._is_frozen = True

        freeze_model(copied)
        return copied

    def unfrozen_copy(self):
        """A copy of this Settings object, which can be modified

        Settings objects become frozen when within a Protocol.  If you *really*
        need to reverse this, this method is how.
        """
        copied = self.model_copy(deep=True)

        def unfreeze_model(model):
            submodels = (
                mod for field in model.model_fields if isinstance(mod := getattr(model, field), SettingsBaseModel)
            )
            for mod in submodels:
                unfreeze_model(mod)

            model._is_frozen = False

        unfreeze_model(copied)

        return copied

    @property
    def is_frozen(self):
        """If this Settings object is frozen and cannot be modified"""
        return self._is_frozen

    def __setattr__(self, name, value):
        if name != "_is_frozen" and self._is_frozen:
            raise AttributeError(
                f"Cannot set '{name}': Settings are immutable once attached"
                " to a Protocol and cannot be modified. Modify Settings "
                "*before* creating the Protocol."
            )
        return super().__setattr__(name, value)

    def __eq__(self, other: Any) -> bool:
        # reproduces pydantic v1 equality, since v2 checks for private attr equality,
        # which results in frozen/unfrozen objects not being equal
        # https://github.com/pydantic/pydantic/blob/2486e068e85c51728c9f2d344cfee2f7e11d555c/pydantic/v1/main.py#L911
        if isinstance(other, _BaseModel):
            return self.model_dump() == other.model_dump()
        else:
            return self.model_dump() == other


class ThermoSettings(SettingsBaseModel):
    """Settings for thermodynamic parameters.

    .. note::
       No checking is done to ensure a valid thermodynamic ensemble is
       possible.
    """

    temperature: KelvinQuantity | None = Field(None, description="Simulation temperature in kelvin)")
    pressure: AtmQuantity | None = Field(None, description="Simulation pressure in standard atmosphere (atm)")
    ph: PositiveFloat | None = Field(None, description="Simulation pH")
    redox_potential: VoltsQuantity | None = Field(None, description="Simulation redox potential in millivolts (mV).")


class BaseForceFieldSettings(SettingsBaseModel, abc.ABC):
    """Base class for ForceFieldSettings objects"""

    ...


def _to_lowercase(value: Any):
    """make any string input lowercase"""
    if isinstance(value, (str)):
        return value.lower()
    else:
        return value


class OpenMMSystemGeneratorFFSettings(BaseForceFieldSettings):
    """Parameters to set up the force field with OpenMM ForceFields

    .. note::
       Currently, this stores what is needed for the
       :class:`openmmforcefields.system_generators.SystemGenerator` signature.
       See the `OpenMMForceField SystemGenerator documentation`_ for more details.


    .. _`OpenMMForceField SystemGenerator documentation`:
       https://github.com/openmm/openmmforcefields#automating-force-field-management-with-systemgenerator
    """

    constraints: Annotated[Literal["hbonds", "allbonds", "hangles"], BeforeValidator(_to_lowercase)] | None = "hbonds"
    """Constraints to be applied to system.
       One of 'hbonds', 'allbonds', 'hangles' or None, default 'hbonds'"""
    rigid_water: bool = True
    """Whether to use a rigid water model. Default True"""
    hydrogen_mass: float = 3.0
    """Mass to be repartitioned to hydrogens from neighbouring
       heavy atoms (in amu), default 3.0"""

    forcefields: list[str] = [
        "amber/ff14SB.xml",  # ff14SB protein force field
        "amber/tip3p_standard.xml",  # TIP3P and recommended monovalent ion parameters
        "amber/tip3p_HFE_multivalent.xml",  # for divalent ions
        "amber/phosaa10.xml",  # Handles THE TPO
    ]
    """List of force field paths for all components except :class:`SmallMoleculeComponent` """

    small_molecule_forcefield: str = "openff-2.1.1"  # other default ideas 'openff-2.0.0', 'gaff-2.11', 'espaloma-0.2.0'
    """Name of the force field to be used for :class:`SmallMoleculeComponent` """

    nonbonded_method: str = "PME"
    """
    Method for treating nonbonded interactions, options are currently
    "CutoffNonPeriodic", "CutoffPeriodic", "Ewald", "LJPME", "NoCutoff", "PME".
    Default PME.
    """
    # TODO: currently, serialization scheme doesn't work for default values since we're not using PlainValidator
    # see https://github.com/pydantic/pydantic/issues/11446
    nonbonded_cutoff: Annotated[NanometerQuantity, Ge(0)] = 1.0 * unit.nanometer  # type: ignore
    """
    Cutoff value for short range nonbonded interactions.
    Default 1.0 * unit.nanometer.
    """

    @field_validator("nonbonded_method", mode="after")
    def allowed_nonbonded_methods(cls, v):
        options = ["CutoffNonPeriodic", "CutoffPeriodic", "Ewald", "LJPME", "NoCutoff", "PME"]
        if v.lower() not in [x.lower() for x in options]:
            errmsg = f"Only {options} are allowed nonbonded_methods"
            raise ValueError(errmsg)
        return v


class Settings(SettingsBaseModel):
    """
    Container for all settings needed by a protocol

    This represents the minimal surface that all settings objects will have.

    Protocols can subclass this to extend this to cater for their additional settings.
    """

    forcefield_settings: InstanceOf[BaseForceFieldSettings]
    thermo_settings: InstanceOf[ThermoSettings]

    @classmethod
    def get_defaults(cls):
        return Settings(
            forcefield_settings=OpenMMSystemGeneratorFFSettings(),
            thermo_settings=ThermoSettings(temperature=300 * unit.kelvin),
        )
