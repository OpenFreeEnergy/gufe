# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from rdkit import Chem
try:
    import openmm.app
except ImportError:
    HAS_OPENMM = False
else:
    HAS_OPENMM = True
from typing import Optional

from gufe import Component
from gufe.custom_typing import RDKitMol, OEMol


class ProteinComponent(Component):
    """Wrapper around a Protein representation.

    This representation is immutable.  If you want to make any modifications,
    do this in an appropriate toolkit then remake this class.
    """
    def __init__(self, *,
                 rdkit: Optional[RDKitMol] = None,
                 openmm_top_and_pos: Optional[tuple] = None,
                 name=""):
        self._rdkit = rdkit
        # yes this is fragile and silly, but it'll do for now
        self._openmm_top, self._openmm_pos = openmm
        self._name = name

    @property
    def name(self):
        return self._name

    @classmethod
    def from_pdbfile(cls, pdbfile: str, name=""):
        raise NotImplementedError()

    @classmethod
    def from_pdbxfile(cls, pdbxfile: str, name=""):
        raise NotImplementedError()

    @classmethod
    def from_rdkit(cls, rdkit: RDKitMol, name=""):
        return cls(rdkit=Chem.Mol(rdkit), name=name)

    def to_rdkit(self) -> RDKitMol:
        return Chem.Mol(self._rdkit)

    @classmethod
    def from_openmm_topology_and_positions(cls, top, pos, name=""):
        """Create from OpenMM topology and positions"""
        return cls(openmm_top_and_pos=(top, pos), name=name)

    def to_openmm_topology_and_positions(self):
        # can't convert from rdkit (presumably what's there) to openmm yet
        if self._openmm_top is None:
            raise AttributeError("OpenMM Topology conversion not possible")
        return self._openmm_top, self._openmm_pos

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
