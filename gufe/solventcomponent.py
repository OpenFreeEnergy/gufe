# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from __future__ import annotations

from openff.units import unit
from typing import Optional, Tuple

from gufe import Component
from gufe.storage.generics import (generic_to_storage_ready,
                                   storage_bytes_to_dict)
from gufe.storage.utils import get_defaults_from_init

_CATIONS = {'Cs', 'K', 'Li', 'Na', 'Rb'}
_ANIONS = {'Cl', 'Br', 'F', 'I'}


# really wanted to make this a dataclass but then can't sort & strip ion input
class SolventComponent(Component):
    """Solvent molecules in a chemical state

    This component represents the abstract idea of the solvent and ions present
    around the other components, rather than a list of specific water molecules
    and their coordinates.  This abstract representation is later made concrete
    by specific MD engine methods.
    """
    _smiles: str
    _positive_ion: Optional[str]
    _negative_ion: Optional[str]
    _neutralize: bool
    _ion_concentration: unit.Quantity

    _storage_path = "setup/components/{md5}.json"

    def __init__(self, *,  # force kwarg usage
                 smiles: str = 'O',
                 positive_ion: Optional[str] = None,
                 negative_ion: Optional[str] = None,
                 neutralize: bool = True,
                 ion_concentration: unit.Quantity = None):
        """
        Parameters
        ----------
        smiles : str, optional
          smiles of the solvent, default 'O' (water)
        positive_ion, negative_ion : str, optional
          the pair of ions which is used to neutralize (if neutralize=True) and
          bring the solvent to the required ionic concentration.  Must be a
          positive and negative monoatomic ions, default `None`
        neutralize : bool, optional
          if the net charge on the chemical state is neutralized by the ions in
          this solvent component.  Default `True`
        ion_concentration : openff-units.unit.Quantity, optional
          ionic concentration required, default `None`
          this must be supplied with units, e.g. "1.5 * unit.molar"

        Examples
        --------
        To create a sodium chloride solution at 0.2 molar concentration::

          >>> s = SolventComponent(position_ion='Na', negative_ion='Cl',
          ...                      ion_concentration=0.2 * unit.molar)

        """
        self._smiles = smiles
        if positive_ion is not None:
            norm = positive_ion.strip('-+').capitalize()
            if norm not in _CATIONS:
                raise ValueError(f"Invalid positive ion, got {positive_ion}")
            positive_ion = norm + '+'
        self._positive_ion = positive_ion
        if negative_ion is not None:
            norm = negative_ion.strip('-+').capitalize()
            if norm not in _ANIONS:
                raise ValueError(f"Invalid negative ion, got {negative_ion}")
            negative_ion = norm + '-'
        self._negative_ion = negative_ion

        self._neutralize = neutralize
        if ion_concentration is not None:
            if (not isinstance(ion_concentration, unit.Quantity) or
               not ion_concentration.is_compatible_with(unit.molar)):
                raise ValueError(f"ion_concentration must be given in units of"
                                 f" concentration, got {ion_concentration}")
            # concentration requires both ions be given
            if ion_concentration > 0:
                if self._negative_ion is None or self._positive_ion is None:
                    raise ValueError("Ions must be given for concentration")
        self._ion_concentration = ion_concentration

    @property
    def name(self) -> str:
        return f"{self.smiles}, {self.positive_ion}, {self.negative_ion}"

    @property
    def smiles(self) -> str:
        """SMILES representation of the solvent molecules"""
        return self._smiles

    @property
    def positive_ion(self) -> Optional[str]:
        """The cation in the solvent state"""
        return self._positive_ion

    @property
    def negative_ion(self) -> Optional[str]:
        """The anion in the solvent state"""
        return self._negative_ion

    @property
    def neutralize(self) -> bool:
        """If the solvent neutralizes the system overall"""
        return self._neutralize

    @property
    def ion_concentration(self) -> unit.Quantity:
        """Concentration of ions in the solvent state"""
        return self._ion_concentration

    @property
    def total_charge(self):
        """Solvents don't have a formal charge defined so this returns None"""
        return None

    def __eq__(self, other):
        if not isinstance(other, SolventComponent):
            return NotImplemented
        return (self.smiles == other.smiles and
                self.positive_ion == other.positive_ion and
                self.negative_ion == other.negative_ion and
                self.ion_concentration == other.ion_concentration)

    def __hash__(self):
        return hash((self.smiles, self.positive_ion, self.negative_ion,
                     self.ion_concentration))

    @classmethod
    def from_dict(cls, d):
        """Deserialize from dict representation"""
        d = d.copy()  # local copy since we modify it
        ion_conc = d['ion_concentration']
        if ion_conc:
            d['ion_concentration'] = unit.Quantity.from_tuple(ion_conc)

        return cls(**d)

    def to_dict(self):
        """For serialization"""
        ion_conc = self.ion_concentration
        if ion_conc:
            ion_conc = ion_conc.to_tuple()

        return {'smiles': self.smiles, 'positive_ion': self.positive_ion,
                'negative_ion': self.negative_ion,
                'ion_concentration': ion_conc,
                'neutralize': self._neutralize}

    def to_storage_ready(self):
        default_keys = get_defaults_from_init(self, self.to_dict())
        return generic_to_storage_ready(self, self._storage_path,
                                        default_keys)

    @classmethod
    def from_storage_bytes(cls, serialized_bytes, load_func):
        dct = storage_bytes_to_dict(serialized_bytes, load_func)
        return cls.from_dict(dct)
