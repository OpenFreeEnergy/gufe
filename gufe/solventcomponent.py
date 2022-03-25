# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import dataclasses
from dataclasses import dataclass
from typing import Tuple


# really wanted to make this a dataclass but then can't sort & strip ion input
class SolventComponent:
    """Solvent molecules

    Used to fully represent the chemical state being expressed.
    """
    def __init__(self, smiles='O', ions=None, concentration=None):
        """
        Parameters
        ----------
        smiles : str, optional
          smiles of the solvent, default 'O' (water)
        ions : list of str, optional
          ions in the system, default None
        concentration : float, optional
          ionic concentration required, default None
        """
        self._smiles = smiles
        if ions is not None:
            # strip and sort so that ('Na', 'Cl-') is equivalent to
            # ('Cl', 'Na+')
            self._ions = tuple(sorted(i.strip('+-') for i in ions))
        else:
            self._ions = tuple()
        self._concentration = concentration

    @property
    def smiles(self):
        return self._smiles

    @property
    def ions(self):
        return self._ions

    @property
    def concentration(self):
        return self._concentration

    def __eq__(self, other):
        # TODO: Class type check?
        return (self.smiles == other.smiles and
                self.ions == other.ions and
                self.concentration == other.concentration)

    def __hash__(self):
        return hash((self.smiles, self.ions, self.concentration))

    def to_dict(self):
        """For serialisation"""
        return {'smiles': self.smiles, 'ions': self.ions,
                'concentration': self.concentration}
