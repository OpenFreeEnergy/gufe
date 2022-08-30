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
            
        mol_topology = openmm_PDBFile.topology
        rd_mol = Mol()
        editable_rdmol = EditableMol(rd_mol)

        # Build Topology
        # Add Atoms
        for atomID, atom in enumerate(mol_topology.topology.atoms()):
            atomID_orig = int(atom.index)
            a = Atom(atom.element.atomic_number)
            a.SetAtomMapNum(atomID_orig)

            a.SetProp("name", atom.name)
            a.SetIntProp("id", atomID_orig)
            a.SetIntProp("formalCharge", 0) # WIP

            a.SetProp("resName", atom.residue.name)
            a.SetIntProp("resId", int(atom.residue.index))

            editable_rdmol.AddAtom(a)


        # Add Bonds
        bondtypes = {1: BondType.SINGLE,
                    2: BondType.DOUBLE,
                    3: BondType.TRIPLE}

        for bond in openmm_PDBFile.topology.bonds():
            bond_order = bond.order # WIP
            editable_rdmol.AddBond(bond.atom1.index, bond.atom2.index)#, bond_order)    

        # Add Additionals
        
        #WIP



        # Set Positions
        rd_mol = editable_rdmol.GetMol()
        positions = list(map(list, openmm_PDBFile.positions._value))
        conf = Conformer(0)
        for atom_id, atom_pos in enumerate(positions):
            conf.SetAtomPosition(atom_id, atom_pos)
        rd_mol.AddConformer(conf)

        # Molecule props
        # Adding nums:
        rd_mol.SetIntProp("NumAtoms", mol_topology.getNumAtoms())
        rd_mol.SetIntProp("NumBonds", mol_topology.getNumBonds())
        rd_mol.SetIntProp("NumChains", mol_topology.getNumChains())

        # dimensions
        pbcVs = list(map(list, mol_topology.getPeriodicBoxVectors()._value))
        #UcellDim = list(map(list, top.getUnitCellDimensions()._value))

        # Sequence Settings
        residue_names = [r.name for r in mol_topology.residues()]
        res_seq = " ".join(residue_names) 
        rd_mol.SetProp("sequence", res_seq)


        # Chains
        rd_mol.SetProp("chain_names", str([c.index for c in mol_topology.chains()]))
        rd_mol.SetProp("chain_resi", str([list([r.index for r in c.residues()]) for c in mol_topology.chains()]))
        # Done
        
        return cls(rdkit=rd_mol, name=name)
    
    @classmethod
    def from_pdbxfile(cls, pdbxfile: str, name=""):
        raise NotImplementedError()

    