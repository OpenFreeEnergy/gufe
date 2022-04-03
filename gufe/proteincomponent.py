# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from rdkit import Chem
try:
    import openmm.app
except ImportError:
    HAS_OPENMM = False
else:
    HAS_OPENMM = True

from gufe import Component
from gufe.custom_typing import RDKitMol, OEMol


class ProteinComponent(Component):
    """Wrapper around a Protein representation.

    This representation is immutable.  If you want to make any modifications,
    do this in an appropriate toolkit then remake this class.
    """
    def __init__(self, rdkit: RDKitMol, name=""):
        self._rdkit = rdkit
        self._openmm_rep = None
        self._name = name

    @property
    def name(self):
        return self._name

    @classmethod
    def from_pdbfile(cls, pdbfile: str, name=""):
        m = Chem.MolFromPDBFile(pdbfile, removeHs=False)
        if m is None:
            raise ValueError(f"RDKit failed to produce a molecule from "
                             "{pdbfile}")
        c = cls(rdkit=m, name=name)
        if HAS_OPENMM:
            # if we can build an openmm representation, do that
            openmm_rep = openmm.app.PDBFile(pdbfile)
            c._openmm_rep = openmm_rep

        return c

    @classmethod
    def from_pdbxfile(cls, pdbxfile: str, name=""):
        raise NotImplementedError()

    @classmethod
    def from_rdkit(cls, rdkit: RDKitMol, name=""):
        return cls(rdkit=Chem.Mol(rdkit), name=name)

    def to_rdkit(self) -> RDKitMol:
        return Chem.Mol(self._rdkit)

    def to_openmm_PDBFile(self):
        return self._openmm_rep

    @classmethod
    def from_openff(cls, offmol, name=""):
        raise NotImplementedError()

    def to_openff(self):
        raise NotImplementedError()

    @classmethod
    def from_openeye(cls, oemol: OEMol, name=""):
        raise NotImplementedError()

    def to_openeye(self) -> OEMol:
        raise NotImplementedError()

    def to_dict(self) -> dict:
        raise NotImplementedError()

    @classmethod
    def from_dict(cls, d: dict):
        raise NotImplementedError()

    def __hash__(self):
        return hash((self.name, Chem.MolToSequence(self._rdkit)))

    def __eq__(self, other):
        return hash(self) == hash(other)

    @property
    def total_charge(self):
        return Chem.GetFormalCharge(self._rdkit)
