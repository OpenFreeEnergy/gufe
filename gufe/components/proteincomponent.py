# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import json, ast
from collections import defaultdict

from openmm import unit as omm_unit

from rdkit import Chem
from rdkit.Chem.rdchem import Mol, Atom, Conformer, EditableMol, BondType

from .explicitmoleculecomponent import ExplicitMoleculeComponent
from .sub_files.pdbfile import PDBFile
from ..molhashing import hashmol, deserialize_numpy, serialize_numpy


bond_types = {  1 : BondType.SINGLE,
                2 : BondType.DOUBLE,
                3 : BondType.TRIPLE ,
               None :  BondType.UNSPECIFIED,
               }

negative_ions = ["CL"]
positive_ions = ["NA", "MG"]       


def assign_correct_prop_type(rd_obj, prop_name, prop_value):
    if(isinstance(prop_value, int)):
        rd_obj.SetIntProp(prop_name, prop_value)
    elif(isinstance(prop_value, float)):
        rd_obj.SetDoubleProp(prop_name, prop_value)
    elif(isinstance(prop_value, bool)):
        rd_obj.SetBoolProp(prop_name, prop_value)
    else:
        rd_obj.SetProp(prop_name, str(prop_value))
        
                            
class ProteinComponent(ExplicitMoleculeComponent):
    """Wrapper around a Protein representation.

    This representation is immutable.  If you want to make any modifications,
    do this in an appropriate toolkit then remake this class.
    """
    
    def to_pdb(self, pdb_file_path:str) -> str:
        """Create a string based on SDF.

        This is the primary serialization mechanism for this class.

        """
        for atom in self._rdkit.GetAtoms():
            #props
            name = str(atom.GetProp("name"))
            resn = str(atom.GetProp("resName"))
            resi = int(atom.GetProp("resId"))
            chainn = str(atom.GetProp("chainName"))
            hetatom = eval(atom.GetProp("hetatom"))
            icode = str(atom.GetProp("insertionCode"))

            # add residue information
            mi = Chem.AtomPDBResidueInfo()
            mi.SetName(name)
            mi.SetResidueName(resn)
            mi.SetResidueNumber(resi)
            mi.SetChainId(chainn)
            mi.SetInsertionCode(icode)
            mi.SetOccupancy(0.0)
            mi.SetTempFactor(0.0)
            mi.SetIsHeteroAtom(hetatom)
            
            atom.SetMonomerInfo(mi)
            
        Chem.MolToPDBFile(self._rdkit, pdb_file_path, flavor=1)
        return pdb_file_path
    
    @classmethod
    def from_pdb_string(cls, pdbfile: str, name=""):
        pass
        
    def to_openmmPDBFile(self):
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
        _residue_atom_map = defaultdict(list)
        histidine_resi_atoms = defaultdict(list)
        _residue_icode = defaultdict(str)
        
        # Add Atoms
        for atom in mol_topology.atoms():
            atomID = int(atom.id)
            atomPosIndex = int(atom.index)
            resn = atom.residue.name
            resi = int(atom.residue.id)
            chainn = str(atom.residue.chain.id)
            chaini = int(atom.residue.chain.index)
            icode = str(atom.residue.insertionCode)
            
            #WIP: get HETATOMS, 
            a = Atom(atom.element.atomic_number)
            a.SetAtomMapNum(atomID)

            a.SetProp("name", atom.name)
            a.SetIntProp("id", atomID)
            a.SetIntProp("_posIndex", atomPosIndex)

            a.SetProp("resName", resn)
            a.SetIntProp("resId", resi)
            a.SetProp("insertionCode", icode)

            a.SetProp("chainName", chainn)
            a.SetIntProp("chainId", chaini)
            
            a.SetProp("hetatom", str(False))
            
            #For histidine fixes
            dict_key = str(resi)+"_"+resn
            if("HIS" ==  atom.residue.name):
                histidine_resi_atoms[dict_key].append(atom.name)
            _residue_atom_map[dict_key].append(atomID)
            if(dict_key not in _residue_icode): _residue_icode[dict_key] = icode
            
            editable_rdmol.AddAtom(a)

        # Add Bonds
        for bond in mol_topology.bonds():
            bond_order = bond_types[bond.order]  
            editable_rdmol.AddBond(beginAtomIdx=bond.atom1.index, endAtomIdx=bond.atom2.index, order=bond_order)    

        # Set Positions
        # WIP: Make multi frame safe
        rd_mol = editable_rdmol.GetMol()
        positions = openmm_PDBFile.positions.value_in_unit(omm_unit.angstrom)
        
        conf = Conformer(0)
        for atom_id, atom_pos in enumerate(positions):
            conf.SetAtomPosition(atom_id, atom_pos) #unit: nm
        rd_mol.AddConformer(conf)


        # Add Additionals
        # Formal Charge
        atoms = rd_mol.GetAtoms()
        netcharge = 0
        for a in atoms:
            atomic_num = a.GetAtomicNum()
            atom_name = a.GetProp("name")
            resn = a.GetProp("resName") 

            connectivity = sum([int(bond.GetBondType()) for bond in a.GetBonds()]) #
            
            default_valence = periodicTable.GetDefaultValence(atomic_num)
            
            # HISTIDINE FIX  resonance
            # Due to the resonance of the Ns in His (which are frequently de/protonating in proteins), there can be bond type changes between ND1-CE1-NE2. 
            if("HIS" == resn and "N" in atom_name and len(atom_name)>1):
                resi = int(a.GetProp("resId"))
                dict_key = str(resi)+"_"+resn

                histidine_atoms = histidine_resi_atoms[dict_key]
                own_prot = atom_name.replace("N", "H") in histidine_atoms
                other_N = list(filter(lambda x: x.startswith("N") and len(x) > 1 and not atom_name== x, histidine_atoms))[0]
                other_prot = other_N.replace("N", "H") in histidine_atoms

                if(own_prot and not other_prot and connectivity != default_valence):
                    #change bond-order
                    bond_change = [bond for bond in a.GetBonds() if("CE1" in (bond.GetBeginAtom().GetProp("name"),
                                                                            bond.GetEndAtom().GetProp("name")))][0]
                    bond_change.SetBondType(bond_types[1])
                    
                    alternate_atom = [atomB for atomB in rd_mol.GetAtoms() if(atomB.GetProp("resId") == str(resi) and atomB.GetProp("name") == str(other_N))][0]
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

        # Molecule props
        # Adding nums:
        rd_mol.SetProp("ofe-name", name)
        rd_mol.SetIntProp("NumAtoms", mol_topology.getNumAtoms())
        rd_mol.SetIntProp("NumBonds", mol_topology.getNumBonds())
        rd_mol.SetIntProp("NumChains", mol_topology.getNumChains())
        rd_mol.SetDoubleProp("NetCharge", netcharge)

        # Chains
        rd_mol.SetProp("chain_names", str([c.id for c in mol_topology.chains()]))
        rd_mol.SetProp("chain_ids", str([c.index for c in mol_topology.chains()]))

        rd_mol.SetProp("_chain_residues", str([[r.index for r in c.residues()] for c in mol_topology.chains()]))

        # Residues
        res_seq = " ".join([r.name for r in mol_topology.residues()])
        rd_mol.SetProp("sequence", res_seq)
        rd_mol.SetProp("_residue_atom_map", str(dict(_residue_atom_map)))
        rd_mol.SetProp("_residue_icode", str(dict(_residue_icode)))

        # Box dimensions
        pbcVs = list(map(list, mol_topology.getPeriodicBoxVectors()._value)) #unit: nm
        unitCellDim = list(map(float, mol_topology.getUnitCellDimensions()._value)) #unit: nm
        rd_mol.SetProp("PeriodicBoxVectors", str(pbcVs))
        rd_mol.SetProp("UnitCellDimensions", str(unitCellDim))

        
        rd_mol.UpdatePropertyCache(strict=True)
        # Done
                
        return cls(rdkit=rd_mol, name=name)
    
    @classmethod
    def from_pdbxfile(cls, pdbxfile: str, name=""):
        raise NotImplementedError()

    @classmethod
    def _from_dict(cls, ser_dict:dict, name=""):
        
        # Mol
        rd_mol = Mol()
        editable_rdmol = EditableMol(rd_mol)
        
        # Add Atoms
        for atom in ser_dict['atoms']:
            atomic_num = int(atom[0])
            atomic_name = atom[1]
            atomic_fc = atom[2]
            atomic_arom = eval(str(atom[3]))
            atomic_ste = atom[4]
            atomic_props = atom[5]

            a = Atom(atomic_num)
            a.SetAtomMapNum(atomic_props["id"])   
            a.SetFormalCharge(atomic_fc)
            a.SetIsAromatic(atomic_arom)
            #a.SetChiralTag(atomic_ste)
            
            for prop_name, prop_value in atomic_props.items():
                assign_correct_prop_type(rd_obj=a, 
                            prop_name=prop_name, 
                            prop_value=prop_value)
                
            editable_rdmol.AddAtom(a)

        # Add Bonds
        for bond in ser_dict['bonds']:
            atomBeginIdx = bond[0]
            atomEndIdx = bond[1]
            bondType = bond[2]
            editable_rdmol.AddBond(beginAtomIdx=atomBeginIdx, 
                                   endAtomIdx=atomEndIdx, 
                                   order=bondType)    
            
        # Set Positions
        rd_mol = editable_rdmol.GetMol()
        positions = ser_dict["conformers"] 
        for i, frame_pos in enumerate(positions):
            frame_pos = deserialize_numpy(frame_pos)
            conf = Conformer(i)
            for atom_id, atom_pos in enumerate(frame_pos):
                conf.SetAtomPosition(atom_id, atom_pos) #unit: nm
            rd_mol.AddConformer(conf)

        #Adding missing bond info
        for bond_id, bond in enumerate(rd_mol.GetBonds()):
            bond_info = ser_dict['bonds'][bond_id]
            bondArom = bond_info[3]
            bondStereo = bond_info[4]
            bondProp = bond_info[5]
            
            bond.SetIsAromatic(bondArom)
            #bond.SetStereo(bondStereo)
            
            for prop_name, prop_value in bondProp.items():
                assign_correct_prop_type(rd_obj=bond, 
                                         prop_name=prop_name, 
                                         prop_value=prop_value)

        # Add Mol Informations
        for mol_prop, mol_value in ser_dict["molecules"].items():
            if(mol_prop == "sequence"):
                mol_value = " ".join(mol_value)    
            assign_correct_prop_type(rd_obj=rd_mol, 
                                        prop_name=mol_prop, 
                                        prop_value=mol_value)
    
        if("name" in ser_dict):
            name = ser_dict["name"]
            
        return cls(rdkit=rd_mol, name=name)        

    def _to_dict(self) -> dict:
               
        # Standards:
        atoms = [
            (atom.GetAtomicNum(),
                atom.GetProp("name"),
                atom.GetFormalCharge(),
                atom.GetIsAromatic(),
                '', 
                atom.GetPropsAsDict()) #stereoCenter
            for atom in self._rdkit.GetAtoms()
        ]
        
        bonds = [
            (bond.GetBeginAtomIdx(), 
             bond.GetEndAtomIdx(),
             bond.GetBondType(),
             bond.GetIsAromatic(), 
             bond.GetStereo() or '',
             bond.GetPropsAsDict())
            for bond in self._rdkit.GetBonds()
        ]
        
        conformers = [
            serialize_numpy(conf.GetPositions()) #.m_as(unit.angstrom)
            for conf in self._rdkit.GetConformers()
        ]
        
        # Additional Information
        #Get Chains
        chains = ast.literal_eval(self._rdkit.GetProp("_chain_residues"))       
        chain_names = ast.literal_eval(self._rdkit.GetProp("chain_names"))
 
        #Residue info
        residue_atom_map = json.loads(self._rdkit.GetProp("_residue_atom_map").replace("'", "\""))
        residue_icode_map = json.loads(self._rdkit.GetProp("_residue_icode").replace("'", "\""))

        residue_name_id = dict([key.split("_") for key in residue_atom_map.keys()])
        residue_sequence = self._rdkit.GetProp("sequence").split()
        
        
        d = {
            'atoms' : atoms,
            'bonds' : bonds,
            'name' : self.name,
            'conformers' : conformers,
            "molecules" : {
                "sequence": residue_sequence,
                "_residue_atom_map": residue_atom_map,
                "_residue_icode": residue_icode_map,
                "residue_name_id": residue_name_id,
                "chain_names": chain_names,
                "_chain_residues": chains
                }
            }

        return d