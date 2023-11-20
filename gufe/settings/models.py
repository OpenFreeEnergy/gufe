# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""
Pydantic models used for storing settings.
"""

import abc
from typing import Optional, Union

from openff.models.models import DefaultModel
from openff.models.types import FloatQuantity
from openff.units import unit

try:
    from pydantic.v1 import (
        Extra,
        Field,
        PositiveFloat,
        validator,
    )
except ImportError:
    from pydantic import (
        Extra,
        Field,
        PositiveFloat,
        validator,
    )


def _find_key_walk(d: dict, target: str, prefix=()):
    # generator for recursively walking Settings for find_field
    for k, v in d.items():
        if k == target:
            yield prefix + (k,)
        elif isinstance(v, dict):
            yield from _find_key_walk(v, target, prefix=prefix + (k,))


class SettingsBaseModel(DefaultModel):
    """Settings and modifications we want for all settings classes."""

    class Config:
        extra = Extra.forbid
        arbitrary_types_allowed = False
        smart_union = True

    def find_field(self, fieldname: str) -> list[tuple[str]]:
        """Return places within this Settings tree where fieldname occurs

        Searches recursive through the settings tree and returns lists of the
        "route" to fields called fieldname.

        E.g. if this class had the field 'foo' both inside the class and also
        within a settings object called 'bar', this function would return::

          [('foo',), ('bar', 'foo')]
        """
        return list(_find_key_walk(self.dict(), fieldname))

    def set_field(self, fieldname, value):
        """Set a parameter to value, searching through the settings tree

        Designed as a convenience, this will search recursively through the
        different Settings objects held by this settings object to find all
        instances of *fieldname*. If only one instance is found, then this is
        set to *value*, otherwise an error is raised as the field being referred
        to is ambiguous.

        Parameters
        ----------
        fieldname : str
          the field you wish to set
        value : Any
          the value to set fieldname to

        Raises
        ------
        ValueError
          if the fieldname appears more than once in the settings tree
        AttributeError
          if the fieldname is never found in the settings tree
        """
        targets = self.find_field(fieldname)

        if not targets:
            raise AttributeError(f"Failed to find '{fieldname}' in Settings")
        if len(targets) > 1:
            # tidy up this state for better error message
            hits = ["'" + '.'.join(vals) + "'" for vals in targets]
            raise ValueError("Ambiguous fieldname, found possibilities at: "
                             f"{', '.join(hits)}")

        target = targets[0]

        roots, attrname = target[:-1], target[-1]
        # e.g. fieldname 'foo' is held at self.bar.baz
        # target = ['bar', 'baz', 'foo']
        # roots = ['bar', 'baz'], attrname = 'foo'
        obj = self
        for r in roots:
            obj = getattr(obj, r)
        setattr(obj, fieldname, value)


class ThermoSettings(SettingsBaseModel):
    """Settings for thermodynamic parameters.

    .. note::
       No checking is done to ensure a valid thermodynamic ensemble is
       possible.
    """

    temperature: FloatQuantity["kelvin"] = Field(
        None, description="Simulation temperature, default units kelvin"
    )
    pressure: FloatQuantity["standard_atmosphere"] = Field(
        None, description="Simulation pressure, default units standard atmosphere (atm)"
    )
    ph: Union[PositiveFloat, None] = Field(None, description="Simulation pH")
    redox_potential: Optional[float] = Field(
        None, description="Simulation redox potential"
    )


class BaseForceFieldSettings(SettingsBaseModel, abc.ABC):
    """Base class for ForceFieldSettings objects"""
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
    constraints: Optional[str] = 'hbonds'
    """Constraints to be applied to system.
       One of 'hbonds', 'allbonds', 'hangles' or None, default 'hbonds'"""

    rigid_water: bool = True
    remove_com: bool = False
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

    @validator('constraints')
    def constraint_check(cls, v):
        allowed = {'hbonds', 'hangles', 'allbonds'}

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
