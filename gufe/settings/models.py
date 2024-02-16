# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""
Pydantic models used for storing settings.
"""

import abc
import pprint
from typing import Optional, Union

from openff.models.models import DefaultModel
from openff.models.types import FloatQuantity
from openff.units import unit

try:
    from pydantic.v1 import Extra, Field, PositiveFloat, PrivateAttr, validator
except ImportError:
    from pydantic import (
        Extra,
        Field,
        PositiveFloat,
        PrivateAttr,
        validator,
    )

import pydantic


class SettingsBaseModel(DefaultModel):
    """Settings and modifications we want for all settings classes."""

    _is_frozen: bool = PrivateAttr(default_factory=lambda: False)

    class Config:
        """
        :noindex:
        """

        extra = pydantic.Extra.forbid
        arbitrary_types_allowed = False
        smart_union = True

    def _ipython_display_(self):
        pprint.pprint(self.dict())

    def frozen_copy(self):
        """A copy of this Settings object which cannot be modified

        This is intended to be used by Protocols to make their stored Settings
        read-only
        """
        copied = self.copy(deep=True)

        def freeze_model(model):
            submodels = (
                mod for field in model.__fields__ if isinstance(mod := getattr(model, field), SettingsBaseModel)
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
        copied = self.copy(deep=True)

        def unfreeze_model(model):
            submodels = (
                mod for field in model.__fields__ if isinstance(mod := getattr(model, field), SettingsBaseModel)
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
                "*before* creating the Protocol.",
            )
        return super().__setattr__(name, value)


class ThermoSettings(SettingsBaseModel):
    """Settings for thermodynamic parameters.

    .. note::
       No checking is done to ensure a valid thermodynamic ensemble is
       possible.
    """

    temperature: FloatQuantity["kelvin"] = Field(None, description="Simulation temperature, default units kelvin")
    pressure: FloatQuantity["standard_atmosphere"] = Field(
        None,
        description="Simulation pressure, default units standard atmosphere (atm)",
    )
    ph: Union[PositiveFloat, None] = Field(None, description="Simulation pH")
    redox_potential: Optional[float] = Field(None, description="Simulation redox potential")


class BaseForceFieldSettings(SettingsBaseModel, abc.ABC):
    """Base class for ForceFieldSettings objects"""

    class Config:
        """:noindex:"""

        pass

    ...


class OpenMMSystemGeneratorFFSettings(BaseForceFieldSettings):
    """Parameters to set up the force field with OpenMM ForceFields

    .. note::
       Right now we just basically just grab what we need for the
       :class:`openmmforcefields.system_generators.SystemGenerator`
       signature. See the `OpenMMForceField SystemGenerator documentation`_
       for more details.


    .. _`OpenMMForceField SystemGenerator documentation`:
       https://github.com/openmm/openmmforcefields#automating-force-field-management-with-systemgenerator
    """

    class Config:
        """:noindex:"""

        pass

    constraints: Optional[str] = "hbonds"
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

    small_molecule_forcefield: str = "openff-2.0.0"  # other default ideas 'openff-2.0.0', 'gaff-2.11', 'espaloma-0.2.0'
    """Name of the force field to be used for :class:`SmallMoleculeComponent` """

    nonbonded_method = "PME"
    """
    Method for treating nonbonded interactions, currently only PME and
    NoCutoff are allowed. Default PME.
    """
    nonbonded_cutoff: FloatQuantity["nanometer"] = 1.0 * unit.nanometer
    """
    Cutoff value for short range nonbonded interactions.
    Default 1.0 * unit.nanometer.
    """

    @validator("nonbonded_method")
    def allowed_nonbonded(cls, v):
        if v.lower() not in ["pme", "nocutoff"]:
            errmsg = "Only PME and NoCutoff are allowed nonbonded_methods"
            raise ValueError(errmsg)
        return v

    @validator("nonbonded_cutoff")
    def is_positive_distance(cls, v):
        # these are time units, not simulation steps
        if not v.is_compatible_with(unit.nanometer):
            raise ValueError("nonbonded_cutoff must be in distance units " "(i.e. nanometers)")
        if v < 0:
            errmsg = "nonbonded_cutoff must be a positive value"
            raise ValueError(errmsg)
        return v

    @validator("constraints")
    def constraint_check(cls, v):
        allowed = {"hbonds", "hangles", "allbonds"}

        if not (v is None or v.lower() in allowed):
            raise ValueError(f"Bad constraints value, use one of {allowed}")

        return v


class Settings(SettingsBaseModel):
    """
    Container for all settings needed by a protocol

    This represents the minimal surface that all settings objects will have.

    Protocols can subclass this to extend this to cater for their additional settings.
    """

    forcefield_settings: BaseForceFieldSettings
    thermo_settings: ThermoSettings

    @classmethod
    def get_defaults(cls):
        return Settings(
            forcefield_settings=OpenMMSystemGeneratorFFSettings(),
            thermo_settings=ThermoSettings(temperature=300 * unit.kelvin),
        )
