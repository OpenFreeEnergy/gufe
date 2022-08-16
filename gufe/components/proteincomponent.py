# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from openff.toolkit.topology import Molecule

from .explicitmoleculecomponent import ExplicitMoleculeComponent
from ..custom_typing import RDKitMol, OEMol



class ProteinComponent(ExplicitMoleculeComponent):
    """Wrapper around a Protein representation.

    This representation is immutable.  If you want to make any modifications,
    do this in an appropriate toolkit then remake this class.
    """
    @classmethod
    def from_pdbfile(cls, pdbfile: str, name=""):
        openff_mol = Molecule.from_polymer_pdb(pdbfile)
        return cls.from_openff(openff_mol, name=name)

    @classmethod
    def from_pdbxfile(cls, pdbxfile: str, name=""):
        raise NotImplementedError()
