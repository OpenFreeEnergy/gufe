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
    def __init__(self, openmm_top, openmm_pos, name=""):
        """
        Parameters
        ----------
        openmm_top : openmm.app.Topology
          the Topology object
        openmm_pos : openmm.unit.Quantity
          the positions for this Topology
        name : str, optional
          identifier for this Protein, used as the hash
        """
        # yes this is fragile and silly, but it'll do for now
        self._openmm_top = openmm_top
        self._openmm_pos = openmm_pos
        self._name = name

    @property
    def name(self):
        return self._name

    @classmethod
    def from_pdbfile(cls, pdbfile: str, name=""):
        if not HAS_OPENMM:
            raise ImportError("OpenMM is currently required")
        f = openmm.app.PDBFile(pdbfile)

        return cls(f.topology, f.positions, name)

    @classmethod
    def from_pdbxfile(cls, pdbxfile: str, name=""):
        raise NotImplementedError()

    @classmethod
    def from_rdkit(cls, rdkit: RDKitMol, name=""):
        raise NotImplementedError()

    def to_rdkit(self) -> RDKitMol:
        raise NotImplementedError()

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

    def _to_dict(self) -> dict:
        raise NotImplementedError()

    @classmethod
    def _from_dict(cls, d: dict):
        raise NotImplementedError()

    @property
    def total_charge(self):
        return None
