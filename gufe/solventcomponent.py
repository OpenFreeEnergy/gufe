# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from openff.units import unit
from typing import Optional, Tuple

from gufe import Component

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
    _ions: Optional[Tuple[str, str]]
    _neutralize: bool
    _ion_concentration: unit.Quantity

    def __init__(self, *,  # force kwarg usage
                 smiles: str = 'O',
                 ions: Tuple[str, str] = None,
                 neutralize: bool = True,
                 ion_concentration: unit.Quantity = None):
        """
        Parameters
        ----------
        smiles : str, optional
          smiles of the solvent, default 'O' (water)
        ions : tuple of str, optional
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

          >>> s = SolventComponent(ions=('Na', 'Cl'),
          ...                      ion_concentration=0.2 * unit.molar)

        """
        self._smiles = smiles
        if ions is not None:
            # normalize: strip, sort and capitalize so that ('Na', 'Cl-') is
            # equivalent to ('Cl', 'Na+')
            norm_ions = tuple(sorted(i.strip('+-').capitalize() for i in ions))
            # check ions make sense
            if len(norm_ions) != 2:
                raise ValueError(f"Must specify exactly two ions, got {ions}")
            n_positive = sum(1 for i in norm_ions if i in _CATIONS)
            n_negative = sum(1 for i in norm_ions if i in _ANIONS)
            if n_positive != 1:
                raise ValueError(f"Must give one positive ion, got {ions}")
            if n_negative != 1:
                raise ValueError(f"Must give one negative ion, got {ions}")
            # mypy gets confused here...
            self._ions = norm_ions  # type: ignore
        else:
            self._ions = None
        self._neutralize = neutralize
        if ion_concentration is not None:
            if (not isinstance(ion_concentration, unit.Quantity) or
               not ion_concentration.is_compatible_with(unit.molar)):
                raise ValueError(f"ion_concentration must be given in units of"
                                 f" concentration, got {ion_concentration}")
        self._ion_concentration = ion_concentration

    @property
    def smiles(self) -> str:
        """SMILES representation of the solvent molecules"""
        return self._smiles

    @property
    def ions(self) -> Optional[Tuple[str, str]]:
        """The ions in the solvent state"""
        return self._ions

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
        try:
            return (self.smiles == other.smiles and
                    self.ions == other.ions and
                    self.ion_concentration == other.ion_concentration)
        except AttributeError:
            return False

    def __hash__(self):
        return hash((self.smiles, self.ions, self.ion_concentration))

    @classmethod
    def from_dict(cls, d):
        """Deserialize from dict representation"""
        return cls(**d)

    def to_dict(self):
        """For serialization"""
        return {'smiles': self.smiles, 'ions': self.ions,
                'ion_concentration': self.ion_concentration,
                'neutralize': self._neutralize}
