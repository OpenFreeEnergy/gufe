# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from typing import Tuple
import math
from openff.units import unit

from gufe import Component


# really wanted to make this a dataclass but then can't sort & strip ion input
class SolventComponent(Component):

    """Solvent molecules in a chemical state

    This component represents the abstract idea of the solvent and ions present
    around the other components, rather than a list of specific water molecules
    and their coordinates.  This abstract representation is later made concrete
    by specific methods.
    """
    def __init__(self, smiles: str = 'O',
                 ions: Tuple[str, ...] = None,
                 neutralize: bool = True,
                 concentration: unit.Quantity = None):
        """
        Parameters
        ----------
        smiles : str, optional
          smiles of the solvent, default 'O' (water)
        ions : list of str, optional
          ions in the system, default `None`
        neutralize : bool, optional
          if the net charge on the chemical state is neutralized by the ions in this
          solvent component.  Default `True`
        concentration : openff-units.unit.Quantity, optional
          ionic concentration required, default `None`
          this must be supplied with units, e.g. "1.5 * unit.molar"
        """
        self._smiles = smiles
        if ions is not None:
            # strip and sort so that ('Na', 'Cl-') is equivalent to
            # ('Cl', 'Na+')
            self._ions = tuple(sorted(i.strip('+-').capitalize()
                                      for i in ions))
        else:
            self._ions = tuple()
        self._neutralize = neutralize
        self._concentration = concentration

    @property
    def smiles(self) -> str:
        """SMILES representation of the solvent molecules"""
        return self._smiles

    @property
    def ions(self) -> Tuple[str, ...]:
        """The ions in the solvent state"""
        return self._ions

    @property
    def neutralize(self) -> bool:
        """If the solvent neutralizes the system overall"""
        return self._neutralize

    @property
    def concentration(self) -> unit.Quantity:
        """Concentration of ions in the solvent state"""
        return self._concentration

    def formal_charge(self) -> int:
        """Solvents don't have a formal charge defined so this returns nan"""
        return math.nan

    def __eq__(self, other):
        try:
            return (self.smiles == other.smiles and
                    self.ions == other.ions and
                    self.concentration == other.concentration)
        except AttributeError:
            return False

    def __hash__(self):
        return hash((self.smiles, self.ions, self.concentration))

    @classmethod
    def from_dict(cls, d):
        """Deserialize from dict representation"""
        return cls(**d)

    def to_dict(self):
        """For serialization"""
        return {'smiles': self.smiles, 'ions': self.ions,
                'concentration': self.concentration,
                'neutralize': self._neutralize}
