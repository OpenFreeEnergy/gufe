# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import ast
import json
import io
import numpy as np
from typing import Union
from collections import defaultdict

from openmm import Vec3
from openmm import app
from openmm import unit as omm_unit
from openmm.app.pdbxfile import PDBxFile

from rdkit import Chem
from rdkit.Chem.rdchem import Mol, Atom, Conformer, EditableMol, BondType

from openff.toolkit.topology import Molecule as OFFMolecule

from .explicitmoleculecomponent import ExplicitMoleculeComponent
from .sub_files.pdbfile import PDBFile
from openmm.app.pdbxfile import PDBxFile

from ..molhashing import deserialize_numpy, serialize_numpy
from ..custom_typing import OEMol


bond_types = {1: BondType.SINGLE,
              2: BondType.DOUBLE,
              3: BondType.TRIPLE,
              None: BondType.UNSPECIFIED,
              }

negative_ions = ["CL"]
positive_ions = ["NA", "MG", 'ZN']


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

        .. note::
        This class is a read-only representation of a protein, if you want to
        edit the molecule do this in an appropriate toolkit **before** creating
        an instance from this class.

        Parameters
        ----------
        rdkit : rdkit.Mol
            rdkit representation of the protein
        name : str, optional
           of the protein, by default ""
    """

    # FROM
    @classmethod
    def from_pdbfile(cls, pdbfile: str, name: str = ""):
        """
        Create ``ProteinComponent`` from PDB-formatted string.

        This is the primary deserialization mechanism for this class.

        Parameters
        ----------
        pdbfile : str
            path to the pdb file.
        name : str, optional
            name of the input protein, by default ""

        Returns
        -------
        ProteinComponent
            the deserialized molecule
        """
        openmm_PDBFile = PDBFile(pdbfile)
        return cls._from_openmmPDBFile(
            openmm_PDBFile=openmm_PDBFile, name=name)

    @classmethod
    def _from_openmmPDBFile(cls, openmm_PDBFile: PDBFile, name: str = ""):
        """
        This Function deserializes openmmPDBFile

        Parameters
        ----------
        openmm_PDBFile : PDBFile
            object of the protein
        name : str
            name of the protein

        Returns
        -------
        ProteinComponent
            the deserialized molecule
        """
        periodicTable = Chem.GetPeriodicTable()
        mol_topology = openmm_PDBFile.getTopology()

        rd_mol = Mol()
        editable_rdmol = EditableMol(rd_mol)

        # Build Topology
        _residue_atom_map = defaultdict(list)
        histidine_resi_atoms = defaultdict(list)
        _residue_icode = defaultdict(str)
        _residue_index = defaultdict(int)
        _residue_id = defaultdict(int)
        # Add Atoms
        for atom in mol_topology.atoms():
            atomID = int(atom.id)
            atomPosIndex = int(atom.index)
            resn = atom.residue.name
            resi = int(atom.residue.id)
            resind = int(atom.residue.index)
            chainn = str(atom.residue.chain.id)
            chaini = int(atom.residue.chain.index)
            icode = str(atom.residue.insertionCode)

            # WIP: get HETATOMS,
            a = Atom(atom.element.atomic_number)
            a.SetAtomMapNum(atomID)

            a.SetProp("name", atom.name)
            a.SetIntProp("id", atomID)
            a.SetIntProp("_posIndex", atomPosIndex)

            a.SetProp("resName", resn)
            a.SetIntProp("resId", resi)
            a.SetIntProp("resInd", resind)
            a.SetProp("insertionCode", icode)

            a.SetProp("chainName", chainn)
            a.SetIntProp("chainId", chaini)

            a.SetProp("hetatom", str(False))

            # For histidine fixes
            dict_key = str(resind) + "_" + resn
            if("HIS" == atom.residue.name):
                histidine_resi_atoms[dict_key].append(atom.name)
            _residue_atom_map[dict_key].append(atomID)
            if(dict_key not in _residue_icode):
                _residue_icode[dict_key] = icode
            if(dict_key not in _residue_index):
                _residue_index[dict_key] = resind
            if(dict_key not in _residue_id):
                _residue_id[dict_key] = resi

            editable_rdmol.AddAtom(a)

        # Add Bonds
        for bond in mol_topology.bonds():
            bond_order = bond_types[bond.order]
            editable_rdmol.AddBond(
                beginAtomIdx=bond.atom1.index,
                endAtomIdx=bond.atom2.index,
                order=bond_order)

        # Set Positions
        rd_mol = editable_rdmol.GetMol()
        positions = np.array(openmm_PDBFile.positions.value_in_unit(
            omm_unit.angstrom), ndmin=3)

        for frame_id, frame in enumerate(positions):
            conf = Conformer(frame_id)
            for atom_id, atom_pos in enumerate(frame):
                conf.SetAtomPosition(atom_id, atom_pos)
            rd_mol.AddConformer(conf)

        # Add Additionals
        # Formal Charge
        atoms = rd_mol.GetAtoms()
        netcharge = 0
        _charged_resi = defaultdict(int)
        for a in atoms:
            atomic_num = a.GetAtomicNum()
            atom_name = a.GetProp("name")
            resn = a.GetProp("resName")
            resind = int(a.GetProp("resInd"))
            dict_key = str(resind) + "_" + resn

            connectivity = sum([int(bond.GetBondType())
                               for bond in a.GetBonds()])

            default_valence = periodicTable.GetDefaultValence(atomic_num)

            # HISTIDINE FIX  resonance
            # Due to the resonance of the Ns in His (which are frequently
            # de/protonating in proteins), there can be bond type changes
            # between ND1-CE1-NE2.
            if("HIS" == resn and "N" in atom_name and atom_name != "N"):
                resind = int(a.GetProp("resInd"))
                dict_key = str(resind) + "_" + resn

                histidine_atoms = histidine_resi_atoms[dict_key]
                own_prot = atom_name.replace("N", "H") in histidine_atoms
                other_N = list(filter(lambda x: x.startswith("N") and len(
                    x) > 1 and not atom_name == x, histidine_atoms))[0]
                other_prot = other_N.replace("N", "H") in histidine_atoms

                if(own_prot and not other_prot and connectivity != default_valence):
                    # change bond-order
                    bond_change = [
                        bond for bond in a.GetBonds() if(
                            "CE1" in (
                                bond.GetBeginAtom().GetProp("name"),
                                bond.GetEndAtom().GetProp("name")))][0]
                    bond_change.SetBondType(bond_types[1])

                    alternate_atom = [atomB for atomB in rd_mol.GetAtoms() if(atomB.GetProp(
                        "resInd") == str(resind) and atomB.GetProp("name") == str(other_N))][0]
                    bond_change = [
                        bond for bond in alternate_atom.GetBonds() if(
                            "CE1" in (
                                bond.GetBeginAtom().GetProp("name"),
                                bond.GetEndAtom().GetProp("name")))][0]
                    bond_change.SetBondType(bond_types[2])

                    # update_new connectivity
                    connectivity = sum([int(bond.GetBondType())
                                        for bond in a.GetBonds()])
            # HISTIDINE FIX DONE

            if(connectivity == 0):  # ions:
                if(atom_name in positive_ions):
                    fc = default_valence  # e.g. Sodium ions
                elif(atom_name in negative_ions):
                    fc = -default_valence  # e.g. Chlorine ions
                else:
                    raise ValueError("I don't know this Ion! \t" + atom_name)
            elif(default_valence > connectivity):
                fc = -(default_valence - connectivity)  # negative charge
            elif(default_valence < connectivity):
                fc = +(connectivity - default_valence)  # positive charge
            else:
                fc = 0  # neutral

            a.SetFormalCharge(fc)
            a.UpdatePropertyCache(strict=True)
            if(fc != 0):
                _charged_resi[dict_key] += fc
            netcharge += fc

        # Molecule props
        # Adding nums:
        rd_mol.SetProp("ofe-name", name)
        rd_mol.SetIntProp("NumAtoms", mol_topology.getNumAtoms())
        rd_mol.SetIntProp("NumBonds", mol_topology.getNumBonds())
        rd_mol.SetIntProp("NumChains", mol_topology.getNumChains())
        rd_mol.SetDoubleProp("NetCharge", netcharge)

        # Chains
        rd_mol.SetProp("chain_names", str(
            [c.id for c in mol_topology.chains()]))
        rd_mol.SetProp("chain_ids", str(
            [c.index for c in mol_topology.chains()]))

        rd_mol.SetProp("_chain_residues", str(
            [[r.index for r in c.residues()] for c in mol_topology.chains()]))

        # Residues
        res_seq = " ".join([r.name for r in mol_topology.residues()])
        rd_mol.SetProp("sequence", res_seq)
        rd_mol.SetProp("_residue_atom_map", str(dict(_residue_atom_map)))
        rd_mol.SetProp("_residue_index", str(dict(_residue_index)))
        rd_mol.SetProp("_residue_id", str(dict(_residue_id)))
        rd_mol.SetProp("_residue_icode", str(dict(_residue_icode)))
        rd_mol.SetProp("_charged_res", str(dict(_charged_resi)))

        # Box dimensions
        pbcVs = mol_topology.getPeriodicBoxVectors()
        if(pbcVs is not None):
            pbcVs = list(map(list, pbcVs.value_in_unit(omm_unit.angstrom)))

        unitCellDim = mol_topology.getUnitCellDimensions()
        if(unitCellDim is not None):
            unitCellDim = list(
                map(float, unitCellDim.value_in_unit(omm_unit.angstrom)))

        rd_mol.SetProp("periodic_box_vectors", str(pbcVs))
        rd_mol.SetProp("unit_cell_dimensions", str(unitCellDim))

        rd_mol.UpdatePropertyCache(strict=True)
        # Done

        return cls(rdkit=rd_mol, name=name)

    @classmethod
    def _from_dict(cls, ser_dict: dict, name: str = ""):
        """Deserialize from dict representation"""

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
            # a.SetChiralTag(atomic_ste)

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
                conf.SetAtomPosition(atom_id, atom_pos)  # unit: nm
            rd_mol.AddConformer(conf)

        # Adding missing bond info
        for bond_id, bond in enumerate(rd_mol.GetBonds()):
            bond_info = ser_dict['bonds'][bond_id]
            bondArom = bond_info[3]
            bondStereo = bond_info[4]
            bondProp = bond_info[5]

            bond.SetIsAromatic(bondArom)
            # bond.SetStereo(bondStereo)

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

    # TO
    def to_openmm_topology(self) -> app.Topology:
        """
        serialize the topology of the protein to openmm.app.Topology

        Returns
        -------
        app.Topology
            resulting topology obj.
        """
        dict_prot = self.to_dict()

        top = app.Topology()

        # Chains
        chains = []
        for chain_name in dict_prot['molecules']["chain_names"]:
            c = top.addChain(id=chain_name)
            chains.append(c)

        # Residues:
        residues = {}
        for res_lab, resind in sorted(
                dict_prot['molecules']["_residue_index"].items(), key=lambda x: x[1]):
            resi, resn = res_lab.split("_")

            resind = dict_prot['molecules']["_residue_index"][res_lab]
            icode = dict_prot['molecules']["_residue_icode"][res_lab]
            resi = int(resi)

            chain_id = int([i for i, v in enumerate(
                dict_prot['molecules']["_chain_residues"]) if(resind in v)][0])
            chain = chains[chain_id]

            # print(resi, resn, chain_id, chain)

            r = top.addResidue(name=resn, id=resind,
                               chain=chain, insertionCode=icode)
            residues.update({chain.id + "_" + str(resi): r})

        # Atoms
        atoms = {}
        for atom in sorted(dict_prot['atoms'], key=lambda x: x[5]["id"]):
            key = atom[5]["chainName"] + "_" + str(atom[5]["resInd"])
            r = residues[key]
            aid = atom[5]["id"]
            atom = top.addAtom(name=atom[1],
                               residue=r,
                               id=aid,
                               element=app.Element.getByAtomicNumber(atom[0])
                               )
            atoms[atom.index] = atom  # true?

        # Bonds
        for bond in dict_prot['bonds']:
            top.addBond(atom1=atoms[bond[0]],
                        atom2=atoms[bond[1]],
                        type=bond[2],
                        order=bond[2])

        # Geometrics
        if(dict_prot['molecules']["unit_cell_dimensions"] != "None"):
            top.setUnitCellDimensions(
                Vec3(*dict_prot['molecules']["unit_cell_dimensions"]) *
                omm_unit.angstrom)

        if(dict_prot['molecules']["periodic_box_vectors"] != "None"):
            top.setPeriodicBoxVectors(
                list(map(lambda x: Vec3(*x), dict_prot['molecules']["periodic_box_vectors"])) *
                omm_unit.angstrom)

        return top

    def to_openmm_positions(self) -> omm_unit.Quantity:
        """
        serialize the positions to openmm.unit.Quantity
        ! only one frame at the moment!

        Returns
        -------
        omm_unit.Quantity
            Quantity containing protein atom positions
        """
        np_pos = deserialize_numpy(self.to_dict()["conformers"][0])
        openmm_pos = list(map(lambda x: Vec3(*x), np_pos)) * omm_unit.angstrom

        return openmm_pos

    def to_pdbFile(self, out_path: Union[str, io.FileIO] = None) -> str:
        """
        serialize protein to pdb file.

        Parameters
        ----------
        out_path : str, optional
            provide path or FileIO to the resulting file, by default None

        Returns
        -------
        str
            string path to the resulting pdb.
        """
        # get top:
        openmm_top = self.to_openmm_topology()

        # get pos:
        openmm_pos = self.to_openmm_positions()

        # write file
        if(isinstance(out_path, str)):
            out_file = open(out_path, "w")
        else:
            out_file = out_path

        PDBFile.writeFile(
            topology=openmm_top,
            positions=openmm_pos,
            file=out_file)

        return out_path

    def to_pdbxFile(self, out_path: str = None) -> str:
        """
            serialize protein to pdbx file.

            Parameters
            ----------
            out_path : str, optional
                provide path or FileIO to the resulting file, by default None

            Returns
            -------
            str
                string path to the resulting pdbx.
        """
        # get top:
        top = self.to_openmm_topology()

        # get pos:
        np_pos = deserialize_numpy(self.to_dict()["conformers"][0])
        openmm_pos = list(map(lambda x: Vec3(*x), np_pos)) * omm_unit.angstrom

        # write file
        if(isinstance(out_path, str)):
            out_file = open(out_path, "w")
        else:
            out_file = out_path

        PDBxFile.writeFile(topology=top, positions=openmm_pos, file=out_file)

        return out_path

    def _to_dict(self) -> dict:
        """Serialize to dict representation"""

        atoms = []
        for atom in self._rdkit.GetAtoms():
            # Standards:
            try:
                name = atom.GetProp("name")
            except KeyError:  # this is default fallback if ff atom name was not stored. mainly used if an rdkit structure is passed.
                name = atom.GetSymbol()

            atoms.append(
                (atom.GetAtomicNum(),
                    name,
                    atom.GetFormalCharge(),
                    atom.GetIsAromatic(),
                    '',
                    atom.GetPropsAsDict()))  # stereoCenter

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
            serialize_numpy(conf.GetPositions())  # .m_as(unit.angstrom)
            for conf in self._rdkit.GetConformers()
        ]

        # Additional Information for the mol:
        molecule_props = {}
        for prop_key, prop_value in self._rdkit.GetPropsAsDict(
                includePrivate=True).items():
            if(prop_key == "sequence"):
                residue_sequence = prop_value.split()
                molecule_props["sequence"] = residue_sequence
            elif(isinstance(prop_value, str) and prop_value.startswith("{")):
                val = json.loads(prop_value.replace("'", "\""))
                molecule_props[prop_key] = val
            elif(isinstance(prop_value, str) and prop_value.startswith("[")):
                val = ast.literal_eval(prop_value)
                molecule_props[prop_key] = val
            else:
                molecule_props[prop_key] = prop_value

        # Result
        d = {
            'atoms': atoms,
            'bonds': bonds,
            'name': self.name,
            'conformers': conformers,
            "molecules": molecule_props
        }

        return d

    # NOT implemented:
    def to_openeye(self) -> OEMol:
        """OEChem representation of this molecule"""
        raise NotImplementedError()

    def to_openff(self):
        raise NotImplementedError()

    def to_gmx(self, out_gmx_top: str, out_gmx_gro: str):
        raise NotImplementedError()

    def to_amber(self, out_path: str):
        raise NotImplementedError()

    @classmethod
    def from_openeye(cls, oemol: OEMol, name: str = ""):
        raise NotImplementedError()

    @classmethod
    def from_openff(cls, openff: OFFMolecule, name: str = ""):
        raise NotImplementedError()

    @classmethod
    def from_gmx(cls, gmx_top: str, gmx_gro: str):
        raise NotImplementedError()

    @classmethod
    def from_amber(cls, name=""):
        raise NotImplementedError()

    @classmethod
    def from_pdbxfile(cls, pdbxfile: str, name=""):
        raise NotImplemented()
        openmm_PDBxFile = PDBxFile(pdbxfile)
        return cls._from_openmmPDBFile(
            openmm_PDBFile=openmm_PDBxFile, name=name)
