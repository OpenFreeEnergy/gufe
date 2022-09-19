# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import ast
import json
import io
import numpy as np
from typing import Union
from collections import defaultdict

from openmm import app
from openmm import unit as omm_unit
from openmm.app.pdbxfile import PDBxFile

from rdkit import Chem
from rdkit.Chem.rdchem import Mol, Atom, Conformer, EditableMol, BondType

from openff.toolkit.topology import Molecule as OFFMolecule

from .explicitmoleculecomponent import ExplicitMoleculeComponent
from ..vendor.pdb_file.pdbfile import PDBFile
from ..vendor.pdb_file.pdbxfile import PDBFile

from ..molhashing import deserialize_numpy, serialize_numpy
from ..custom_typing import OEMol


bond_types = {
    1: BondType.SINGLE,
    2: BondType.DOUBLE,
    3: BondType.TRIPLE,
    None: BondType.UNSPECIFIED,
}

negative_ions = ["F", "CL", "Br", "I"]
positive_ions = ["NA", "MG", "ZN"]


def assign_correct_prop_type(rd_obj, prop_name, prop_value):
    if isinstance(prop_value, int):
        rd_obj.SetIntProp(prop_name, prop_value)
    elif isinstance(prop_value, float):
        rd_obj.SetDoubleProp(prop_name, prop_value)
    elif isinstance(prop_value, bool):
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
            openmm_PDBFile=openmm_PDBFile, name=name
        )

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
            ishetatom = False

            # WIP: get HETATOMS,
            a = Atom(atom.element.atomic_number)
            a.SetAtomMapNum(atomID)

            a.SetIntProp("id", atomID)
            a.SetIntProp("resId", resi)
            a.SetIntProp("_posIndex", atomPosIndex)

            atom_monomerInfo = Chem.AtomPDBResidueInfo()
            atom_monomerInfo.SetChainId(chainn)
            atom_monomerInfo.SetSegmentNumber(chaini)
            atom_monomerInfo.SetInsertionCode(icode)
            atom_monomerInfo.SetName(atom.name)
            atom_monomerInfo.SetResidueName(resn)
            atom_monomerInfo.SetResidueNumber(resind)
            atom_monomerInfo.SetIsHeteroAtom(ishetatom)

            a.SetMonomerInfo(atom_monomerInfo)

            # additonally possible:
            # mi.SetSerialNumber
            # mi.SetSecondaryStructure
            # mi.SetMonomerType
            # mi.SetAltLoc

            # For molecule props
            dict_key = str(resind) + "_" + resn
            _residue_atom_map[dict_key].append(atomID)
            if dict_key not in _residue_icode:
                _residue_icode[dict_key] = icode
            if dict_key not in _residue_index:
                _residue_index[dict_key] = resind
            if dict_key not in _residue_id:
                _residue_id[dict_key] = resi

            editable_rdmol.AddAtom(a)

        # Add Bonds
        for bond in mol_topology.bonds():
            bond_order = bond_types[bond.order]
            editable_rdmol.AddBond(
                beginAtomIdx=bond.atom1.index,
                endAtomIdx=bond.atom2.index,
                order=bond_order,
            )

        # Set Positions
        rd_mol = editable_rdmol.GetMol()
        positions = np.array(
            openmm_PDBFile.positions.value_in_unit(omm_unit.angstrom), ndmin=3
        )

        for frame_id, frame in enumerate(positions):
            conf = Conformer(frame_id)
            for atom_id, atom_pos in enumerate(frame):
                conf.SetAtomPosition(atom_id, atom_pos)
            rd_mol.AddConformer(conf)

        # Add Additionals
        # Formal Charge
        atoms = rd_mol.GetAtoms()
        netcharge = 0
        _charged_resi: defaultdict = defaultdict(int)
        for a in atoms:
            atomic_num = a.GetAtomicNum()
            atom_name = a.GetMonomerInfo().GetName()
            resn = a.GetMonomerInfo().GetResidueName()
            resind = int(a.GetMonomerInfo().GetResidueNumber())
            dict_key = str(resind) + "_" + resn

            connectivity = sum(
                [int(bond.GetBondType()) for bond in a.GetBonds()]
            )
            default_valence = periodicTable.GetDefaultValence(atomic_num)

            if connectivity == 0:  # ions:
                if atom_name in positive_ions:
                    fc = default_valence  # e.g. Sodium ions
                elif atom_name in negative_ions:
                    fc = -default_valence  # e.g. Chlorine ions
                else:
                    raise ValueError("I don't know this Ion! \t" + atom_name)
            elif default_valence > connectivity:
                fc = -(default_valence - connectivity)  # negative charge
            elif default_valence < connectivity:
                fc = +(connectivity - default_valence)  # positive charge
            else:
                fc = 0  # neutral

            a.SetFormalCharge(fc)
            a.UpdatePropertyCache(strict=True)
            if fc != 0:
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
        rd_mol.SetProp(
            "chain_names", str([c.id for c in mol_topology.chains()])
        )
        rd_mol.SetProp(
            "chain_ids", str([c.index for c in mol_topology.chains()])
        )

        rd_mol.SetProp(
            "_chain_residues",
            str(
                [
                    [r.index for r in c.residues()]
                    for c in mol_topology.chains()
                ]
            ),
        )

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
        if pbcVs is not None:
            pbcVs = list(map(list, pbcVs.value_in_unit(omm_unit.angstrom)))

        unitCellDim = mol_topology.getUnitCellDimensions()
        if unitCellDim is not None:
            unitCellDim = list(
                map(float, unitCellDim.value_in_unit(omm_unit.angstrom))
            )

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
        for atom in ser_dict["atoms"]:
            atomic_num = int(atom[0])
            atomic_name = atom[1]
            atomic_fc = atom[2]
            atomic_arom = eval(str(atom[3]))
            atomic_ste = atom[4]
            atomic_props = atom[5]
            atom_mi_dict = atom[6]

            a = Atom(atomic_num)
            a.SetAtomMapNum(atomic_props["id"])
            a.SetFormalCharge(atomic_fc)
            a.SetIsAromatic(atomic_arom)
            # a.SetChiralTag(atomic_ste)

            # put mi_dict back to class
            atom_monomerInfo = Chem.AtomPDBResidueInfo()
            for key, val in atom_mi_dict.items():
                f = getattr(atom_monomerInfo, "Set" + str(key))
                f(val)

            a.SetMonomerInfo(atom_monomerInfo)

            for prop_name, prop_value in atomic_props.items():
                assign_correct_prop_type(
                    rd_obj=a, prop_name=prop_name, prop_value=prop_value
                )

            editable_rdmol.AddAtom(a)

        # Add Bonds
        for bond in ser_dict["bonds"]:
            atomBeginIdx = bond[0]
            atomEndIdx = bond[1]
            bondType = bond[2]
            editable_rdmol.AddBond(
                beginAtomIdx=atomBeginIdx,
                endAtomIdx=atomEndIdx,
                order=bondType,
            )

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
            bond_info = ser_dict["bonds"][bond_id]
            bondArom = bond_info[3]
            bondStereo = bond_info[4]
            bondProp = bond_info[5]

            bond.SetIsAromatic(bondArom)
            # bond.SetStereo(bondStereo)

            for prop_name, prop_value in bondProp.items():
                assign_correct_prop_type(
                    rd_obj=bond, prop_name=prop_name, prop_value=prop_value
                )

        # Add Mol Informations
        for mol_prop, mol_value in ser_dict["molecules"].items():
            if mol_prop == "sequence":
                mol_value = " ".join(mol_value)
            assign_correct_prop_type(
                rd_obj=rd_mol, prop_name=mol_prop, prop_value=mol_value
            )

        if "name" in ser_dict:
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
        for chain_name in dict_prot["molecules"]["chain_names"]:
            c = top.addChain(id=chain_name)
            chains.append(c)

        # Residues:
        residues = {}
        for res_lab, resind in sorted(
            dict_prot["molecules"]["_residue_index"].items(),
            key=lambda x: x[1],
        ):
            resi, resn = res_lab.split("_")

            resind = dict_prot["molecules"]["_residue_index"][res_lab]
            icode = dict_prot["molecules"]["_residue_icode"][res_lab]
            resi = int(resi)

            chain_id = int(
                [
                    i
                    for i, v in enumerate(
                        dict_prot["molecules"]["_chain_residues"]
                    )
                    if (resind in v)
                ][0]
            )
            chain = chains[chain_id]

            # print(resi, resn, chain_id, chain)

            r = top.addResidue(
                name=resn, id=resind, chain=chain, insertionCode=icode
            )
            residues.update({chain.id + "_" + str(resi): r})

        # Atoms
        atoms = {}
        for atom in sorted(dict_prot["atoms"], key=lambda x: x[5]["id"]):
            aid = atom[5]["id"]
            atom_mi_dict = atom[6]
            chainn = atom_mi_dict["ChainId"]
            resInd = atom_mi_dict["ResidueNumber"]

            key = str(chainn) + "_" + str(resInd)
            r = residues[key]

            atom = top.addAtom(
                name=atom[1],
                residue=r,
                id=aid,
                element=app.Element.getByAtomicNumber(atom[0]),
            )
            atoms[atom.index] = atom

        # Bonds
        for bond in dict_prot["bonds"]:
            top.addBond(
                atom1=atoms[bond[0]],
                atom2=atoms[bond[1]],
                type=bond[2],
                order=bond[2],
            )

        # Geometrics
        if dict_prot["molecules"]["unit_cell_dimensions"] != "None":
            top.setUnitCellDimensions(
                np.array(dict_prot["molecules"]["unit_cell_dimensions"])
                * omm_unit.angstrom
            )
        else:
            top.setUnitCellDimensions(None)

        if dict_prot["molecules"]["periodic_box_vectors"] != "None":
            top.setPeriodicBoxVectors(
                list(
                    map(
                        lambda x: np.array(x),
                        dict_prot["molecules"]["periodic_box_vectors"],
                    )
                )
                * omm_unit.angstrom
            )
        else:
            top.setPeriodicBoxVectors(None)

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
        openmm_pos = (
            list(map(lambda x: np.array(x), np_pos)) * omm_unit.angstrom
        )

        return openmm_pos

    def to_pdbFile(self, out_path: Union[str, io.TextIOBase] = None) -> str:
        """
        serialize protein to pdb file.

        Parameters
        ----------
        out_path : Union[str, io.TextIOBase]
            provide path or any string based stream (e.g. FileIO ) to the resulting file, by default None

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
        if isinstance(out_path, str):
            out_file = open(out_path, "w")
        elif isinstance(out_path, io.TextIOBase):
            out_file = out_path
            out_path = str(out_file.name)
        else:
            raise ValueError("Out path type was not as expected!")

        PDBFile.writeFile(
            topology=openmm_top, positions=openmm_pos, file=out_file
        )

        return out_path

    def to_pdbxFile(
        self, out_path: Union[str, io.TextIOBase] = None
    ) -> str:
        """
        serialize protein to pdbx file.

        Parameters
        ----------
        out_path : str
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
        openmm_pos = (
            list(map(lambda x: np.array(x), np_pos)) * omm_unit.angstrom
        )

        # write file
        if isinstance(out_path, str):
            out_file = open(out_path, "w")
        elif isinstance(out_path, io.TextIOBase):
            out_file = out_path
            out_path = str(out_file.name)
        else:
            raise ValueError("Out path type was not as expected!")

        PDBxFile.writeFile(topology=top, positions=openmm_pos, file=out_file)

        return out_path

    def _to_dict(self) -> dict:
        """Serialize to dict representation"""

        atoms = []
        for atom in self._rdkit.GetAtoms():
            # Standards:
            name = atom.GetMonomerInfo().GetName()

            if name == "":
                name = atom.GetSymbol()

            # collapse monomer info to dict
            atom_monomer_info = atom.GetMonomerInfo()
            mi_dict = {}
            for f in dir(atom_monomer_info):
                if "Get" in f:
                    val = getattr(atom_monomer_info, f)()
                    mi_dict[f.replace("Get", "")] = val

            atoms.append(
                (
                    atom.GetAtomicNum(),
                    name,
                    atom.GetFormalCharge(),
                    atom.GetIsAromatic(),
                    "",  # Stereocent in Smallcomponent
                    atom.GetPropsAsDict(),
                    mi_dict,
                )
            )

        bonds = [
            (
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                bond.GetBondType(),
                bond.GetIsAromatic(),
                bond.GetStereo() or "",
                bond.GetPropsAsDict(),
            )
            for bond in self._rdkit.GetBonds()
        ]

        conformers = [
            serialize_numpy(conf.GetPositions())  # .m_as(unit.angstrom)
            for conf in self._rdkit.GetConformers()
        ]

        # Additional Information for the mol:
        molecule_props = {}
        for prop_key, prop_value in self._rdkit.GetPropsAsDict(
            includePrivate=True
        ).items():
            if prop_key == "sequence":
                residue_sequence = prop_value.split()
                molecule_props["sequence"] = residue_sequence
            elif isinstance(prop_value, str) and prop_value.startswith("{"):
                val = json.loads(prop_value.replace("'", '"'))
                molecule_props[prop_key] = val
            elif isinstance(prop_value, str) and prop_value.startswith("["):
                val = ast.literal_eval(prop_value)
                molecule_props[prop_key] = val
            else:
                molecule_props[prop_key] = prop_value

        # Result
        d = {
            "atoms": atoms,
            "bonds": bonds,
            "name": self.name,
            "conformers": conformers,
            "molecules": molecule_props,
        }

        return d

    # NOT implemented:
    def to_openeye(self) -> OEMol:
        """OEChem representation of this molecule"""
        raise NotImplementedError

    def to_openff(self):
        raise NotImplementedError

    def to_gmx(self, out_gmx_top: str, out_gmx_gro: str):
        raise NotImplementedError

    def to_amber(self, out_path: str):
        raise NotImplementedError

    @classmethod
    def from_openeye(cls, oemol: OEMol, name: str = ""):
        raise NotImplementedError

    @classmethod
    def from_openff(cls, openff: OFFMolecule, name: str = ""):
        raise NotImplementedError

    @classmethod
    def from_gmx(cls, gmx_top: str, gmx_gro: str):
        raise NotImplementedError

    @classmethod
    def from_amber(cls, name=""):
        raise NotImplementedError

    @classmethod
    def from_pdbxfile(cls, pdbxfile: str, name=""):
        raise NotImplemented
        # openmm_PDBxFile = PDBxFile(pdbxfile)
        # return cls._from_openmmPDBFile(
        #    openmm_PDBFile=openmm_PDBxFile, name=name)
