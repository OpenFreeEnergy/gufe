# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from openff.toolkit.topology import Molecule

from rdkit.Chem.rdchem import Mol, Atom, Conformer, EditableMol
from openmm.app import PDBFile

from gufe.components.explicitmoleculecomponent import ExplicitMoleculeComponent



class ProteinComponent(ExplicitMoleculeComponent):
    """Wrapper around a Protein representation.

    This representation is immutable.  If you want to make any modifications,
    do this in an appropriate toolkit then remake this class.
    """
    
    def to_pdb(self) -> str:
        """Create a string based on SDF.

        This is the primary serialization mechanism for this class.

        """
        pass
    
    @classmethod
    def from_pdb_string(cls, pdbfile: str, name=""):
        pass
        
    
    @classmethod
    def from_pdbfile(cls, pdbfile: str, name=""):
        """
        _summary_

        Parameters
        ----------
        pdbfile : str
            _description_
        name : str, optional
            _description_, by default ""

        Returns
        -------
        _type_
            _description_
        """
        openmm_PDBFile = PDBFile(pdbfile)
        return cls.from_openmmPDBFile(openmm_PDBFile=openmm_PDBFile, name=name)
    
    @classmethod
    def from_openmmPDBFile(cls, openmm_PDBFile:PDBFile, name:str):
        """
        This Function serializes openmmPDBFile

        Parameters
        ----------
        openmm_PDBFile : PDBFile
            _description_
        name : str
            _description_

        Returns
        -------
        _type_
            _description_
        """
        rdmol = Mol()
        rdemol = EditableMol(rdmol)

        # Build Topology
        # Add Atoms
        for atom in openmm_PDBFile.topology.atoms():
            a = Atom(atom.element.atomic_number)
            rdemol.AddAtom(a)

        # Add Bonds


        # Add Additionals


        # Set Positions
        rd_mol = rdemol.GetMol()
        positions = list(map(list, openmm_PDBFile.positions._value))
        conf = Conformer(0)

        for atom_id, atom_pos in enumerate(positions):
            conf.SetAtomPosition(atom_id, atom_pos)
            
        rd_mol.AddConformer(conf)

        # Done
        return cls(rd_mol, name=name)
    
    @classmethod
    def from_pdbxfile(cls, pdbxfile: str, name=""):
        raise NotImplementedError()

    