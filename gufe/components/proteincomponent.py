# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from openff.toolkit.topology import Molecule

from rdkit import Chem
from rdkit.Chem.rdchem import Mol, Atom, Conformer, EditableMol, BondType
from openmm.app import PDBFile

from gufe.components.explicitmoleculecomponent import ExplicitMoleculeComponent
from .subfiles.PDBFile import PDBFile


bond_types = {  1 : BondType.SINGLE,
                2 : BondType.DOUBLE,
                3 : BondType.TRIPLE ,
               None :  BondType.UNSPECIFIED,
               }
          
                    
class ProteinComponent(ExplicitMoleculeComponent):
    """Wrapper around a Protein representation.

    This representation is immutable.  If you want to make any modifications,
    do this in an appropriate toolkit then remake this class.
    """
    
    def to_pdb(self, pdb_file_path:str) -> str:
        """Create a string based on SDF.

        This is the primary serialization mechanism for this class.

        """
        Chem.MolToPDBFile(pdb_file_path)
        return pdb_file_path
    
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
        return cls._from_openmmPDBFile(openmm_PDBFile=openmm_PDBFile, name=name)
    
    
    @classmethod
    def _from_openmmPDBFile(cls, openmm_PDBFile:PDBFile, name:str):
        """
        This Function serializes openmmPDBFile to 
        AA - Protonations
        
        Test:        
         - 1.5 serialization test
         - check out files
         - check obj

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
            
                    
        periodicTable = Chem.GetPeriodicTable()
        mol_topology = openmm_PDBFile.topology

        rd_mol = Mol()
        editable_rdmol = EditableMol(rd_mol)

        # Build Topology
        # Add Atoms
        for atom in mol_topology.atoms():
            atomID_orig = int(atom.index)
            a = Atom(atom.element.atomic_number)
            a.SetAtomMapNum(atomID_orig)

            a.SetProp("name", atom.name)
            a.SetIntProp("id", atomID_orig)

            a.SetProp("resName", atom.residue.name)
            a.SetIntProp("resId", int(atom.residue.index))

            editable_rdmol.AddAtom(a)

        # Add Bonds
        for bond in mol_topology.bonds():
            bond_order = bond_types[bond.order]
            editable_rdmol.AddBond(beginAtomIdx=bond.atom1.index, endAtomIdx=bond.atom2.index,)# order=bond_order)    

        # Set Positions
        rd_mol = editable_rdmol.GetMol()
        positions = list(map(list, openmm_PDBFile.positions._value))
        conf = Conformer(0)
        for atom_id, atom_pos in enumerate(positions):
            conf.SetAtomPosition(atom_id, atom_pos) #unit: nm
        rd_mol.AddConformer(conf)

        # Molecule props
        # Adding nums:
        rd_mol.SetIntProp("NumAtoms", mol_topology.getNumAtoms())
        rd_mol.SetIntProp("NumBonds", mol_topology.getNumBonds())
        rd_mol.SetIntProp("NumChains", mol_topology.getNumChains())

        # dimensions
        pbcVs = list(map(list, mol_topology.getPeriodicBoxVectors()._value)) #unit: nm
        unitCellDim = list(map(float, mol_topology.getUnitCellDimensions()._value)) #unit: nm
        rd_mol.SetProp("PeriodicBoxVectors", str(pbcVs))
        rd_mol.SetProp("UnitCellDimensions", str(unitCellDim))

        # Sequence Settings
        residue_names = [r.name for r in mol_topology.residues()]
        res_seq = " ".join(residue_names) 
        rd_mol.SetProp("sequence", res_seq)

        # Chains
        rd_mol.SetProp("chain_names", str([c.index for c in mol_topology.chains()]))
        rd_mol.SetProp("chain_resi", str([list([r.index for r in c.residues()]) for c in mol_topology.chains()]))


        # Add Additionals
        # Set Bondorder
        # WIP: if not possible above
        
        # Formal Charge
        # WIP: I need Bondorder here!
        atoms = rd_mol.GetAtoms()
        netcharge = 0
        for a in atoms:
            a.UpdatePropertyCache()
            connectivity = len(a.GetBonds()) #a.GetTotalValence()
            atomic_num = a.GetAtomicNum()
            default_valence = periodicTable.GetDefaultValence(atomic_num)
            if(default_valence > connectivity):
                fc = -(default_valence-connectivity) # negative charge
            elif(default_valence < connectivity):
                fc = +(connectivity-default_valence) # positive charge
            else:
                fc = 0 # neutral
            
            a.SetFormalCharge(fc)
            #print(connectivity, default_valence, periodicTable.GetElementSymbol(atomic_num))

            #if(fc > 0):
            #    print(periodicTable.GetElementSymbol(atomic_num), fc, connectivity)
            netcharge+=fc
            
        rd_mol.SetDoubleProp("NetCharge", netcharge)
        rd_mol.UpdatePropertyCache()
        # Done
        
        return cls(rdkit=rd_mol, name=name)
    
    @classmethod
    def from_pdbxfile(cls, pdbxfile: str, name=""):
        raise NotImplementedError()

    