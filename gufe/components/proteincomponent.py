# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from openff.toolkit.topology import Molecule

from rdkit import Chem
from rdkit.Chem.rdchem import Mol, Atom, Conformer, EditableMol, BondType
from openff.units import unit

from gufe.components.explicitmoleculecomponent import ExplicitMoleculeComponent
from .sub_files.pdbfile import PDBFile


bond_types = {  1 : BondType.SINGLE,
                2 : BondType.DOUBLE,
                3 : BondType.TRIPLE ,
               None :  BondType.UNSPECIFIED,
               }

negative_ions = ["CL"]
positive_ions = ["NA", "MG"]       
                    
class ProteinComponent(ExplicitMoleculeComponent):
    """Wrapper around a Protein representation.

    This representation is immutable.  If you want to make any modifications,
    do this in an appropriate toolkit then remake this class.
    """
    
    def to_pdb(self, pdb_file_path:str) -> str:
        """Create a string based on SDF.

        This is the primary serialization mechanism for this class.

        """
        Chem.MolToPDBFile(self._rdkit, pdb_file_path)
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
        histidine_resi_atoms = {}
        # Add Atoms
        for atom in mol_topology.atoms():
            atomID_orig = int(atom.index)
            a = Atom(atom.element.atomic_number)
            a.SetAtomMapNum(atomID_orig)

            a.SetProp("name", atom.name)
            a.SetIntProp("id", atomID_orig)

            a.SetProp("resName", atom.residue.name)
            a.SetIntProp("resId", int(atom.residue.index))

            if("HIS" ==  atom.residue.name):
                if(int(atom.residue.index) in histidine_resi_atoms):
                    histidine_resi_atoms[int(atom.residue.index)].append(atom.name)
                else:
                    histidine_resi_atoms[int(atom.residue.index)] = [atom.name]
                
            editable_rdmol.AddAtom(a)

        # Add Bonds
        for bond in mol_topology.bonds():
            bond_order = bond_types[bond.order]  
            editable_rdmol.AddBond(beginAtomIdx=bond.atom1.index, endAtomIdx=bond.atom2.index, order=bond_order)    

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
        rd_mol.SetProp("chain_resi", str([[r.index for r in c.residues()] for c in mol_topology.chains()]))


        # Add Additionals
        # Formal Charge
        atoms = rd_mol.GetAtoms()
        netcharge = 0
        for a in atoms:
            atomic_num = a.GetAtomicNum()
            atom_name = a.GetProp("name")

            connectivity = sum([int(bond.GetBondType()) for bond in a.GetBonds()]) #
            
            default_valence = periodicTable.GetDefaultValence(atomic_num)
            
            # HISTIDINE FIX  resonance
            # Due to the resonance of the Ns in His (which are frequently de/protonating in proteins), there can be bond type changes between ND1-CE1-NE2. 
            if("HIS" == a.GetProp("resName") and "N" in a.GetProp("name") and len(a.GetProp("name"))>1):
                resi = int(a.GetProp("resId"))

                histidine_atoms = histidine_resi_atoms[resi]
                own_prot = a.GetProp("name").replace("N", "H") in histidine_atoms
                other_N = list(filter(lambda x: x.startswith("N") and len(x) > 1 and not atom_name== x, histidine_atoms))[0]
                other_prot = other_N.replace("N", "H") in histidine_atoms

                if(own_prot and not other_prot and connectivity != default_valence):
                    #change bond-order
                    bond_change = [bond for bond in a.GetBonds() if("CE1" in (bond.GetBeginAtom().GetProp("name"),
                                                                            bond.GetEndAtom().GetProp("name")))][0]
                    bond_change.SetBondType(bond_types[1])
                    
                    alternate_atom = [a for a in rd_mol.GetAtoms() if(a.GetProp("resId") == str(resi) and a.GetProp("name") == str(other_N))][0]
                    bond_change = [bond for bond in alternate_atom.GetBonds() if("CE1" in (bond.GetBeginAtom().GetProp("name"),
                                                                                        bond.GetEndAtom().GetProp("name")))][0]
                    bond_change.SetBondType(bond_types[2])  
                connectivity = sum([int(bond.GetBondType()) for bond in a.GetBonds()])

            ### HISTIDINE FIX DONE
            
            if(connectivity == 0): #ions:
                if(atom_name in positive_ions):
                    fc = default_valence  #e.g. Sodium ions
                elif(atom_name in negative_ions):
                    fc = -default_valence  #e.g. Chlorine ions
                else:
                    raise ValueError("I don't know this Ion! \t"+atom_name)     
            elif(default_valence > connectivity):
                fc = -(default_valence-connectivity) # negative charge
            elif(default_valence < connectivity):
                fc = +(connectivity-default_valence) # positive charge
            else:
                fc = 0 # neutral

            a.SetFormalCharge(fc)
            a.UpdatePropertyCache(strict=True)
            
            netcharge+=fc
            
        rd_mol.SetDoubleProp("NetCharge", netcharge)
        rd_mol.UpdatePropertyCache(strict=True)
        # Done
        
        return cls(rdkit=rd_mol, name=name)
    
    @classmethod
    def from_pdbxfile(cls, pdbxfile: str, name=""):
        raise NotImplementedError()

    