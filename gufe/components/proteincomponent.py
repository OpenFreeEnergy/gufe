# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import ast
import io
import json
import string
from collections import defaultdict
from os import PathLike
from typing import Optional, Union

import numpy as np
from openmm import app
from openmm import unit as omm_unit
from rdkit import Chem
from rdkit.Chem.rdchem import Atom, BondType, Conformer, EditableMol, Mol

from ..custom_typing import RDKitMol
from ..molhashing import deserialize_numpy, serialize_numpy
from ..vendor.pdb_file.pdbfile import PDBFile
from ..vendor.pdb_file.pdbxfile import PDBxFile
from .explicitmoleculecomponent import ExplicitMoleculeComponent

_BONDORDERS_OPENMM_TO_RDKIT = {
    1: BondType.SINGLE,
    2: BondType.DOUBLE,
    3: BondType.TRIPLE,
    None: BondType.UNSPECIFIED,
}
_BONDTYPES_OPENMM_TO_RDKIT = {
    app.Single: BondType.SINGLE,
    app.Double: BondType.DOUBLE,
    app.Triple: BondType.TRIPLE,
    app.Aromatic: BondType.AROMATIC,
    None: BondType.UNSPECIFIED,
}
_BONDORDERS_RDKIT_TO_OPENMM = {v: k for k, v in _BONDORDERS_OPENMM_TO_RDKIT.items()}
_BONDTYPES_RDKIT_TO_OPENMM = {v: k for k, v in _BONDTYPES_OPENMM_TO_RDKIT.items()}
_BONDORDER_TO_ORDER = {
    BondType.UNSPECIFIED: 1,  # assumption
    BondType.SINGLE: 1,
    BondType.DOUBLE: 2,
    BondType.TRIPLE: 3,
}


# builtin dict of strings to enum members, boy I hope this is stable
_BONDORDER_STR_TO_RDKIT = Chem.BondType.names
_BONDORDER_RDKIT_TO_STR = {v: k for k, v in _BONDORDER_STR_TO_RDKIT.items()}

_CHIRALITY_RDKIT_TO_STR = {
    Chem.CHI_TETRAHEDRAL_CW: "CW",
    Chem.CHI_TETRAHEDRAL_CCW: "CCW",
    Chem.CHI_UNSPECIFIED: "U",
}
_CHIRALITY_STR_TO_RDKIT = {v: k for k, v in _CHIRALITY_RDKIT_TO_STR.items()}


ions_dict = {
    # Alkali metals
    "LI": 1,
    "NA": 1,
    "K": 1,
    "RB": 1,
    "CS": 1,
    "K+": 1,
    "Na+": 1,
    # Alkaline earth metals
    "Be": 2,
    "MG": 2,
    "CA": 2,
    "SR": 2,
    "BA": 2,
    "Ra": 2,
    # Transition metals
    "CE": 3,
    "Ce": 4,
    "CR": 3,
    "Cr": 2,
    "MN": 2,
    "FE": 3,
    "FE2": 2,
    "CO": 2,
    "NI": 2,
    "CU": 2,
    "CU1": 1,
    "ZN": 2,
    "AG": 1,
    "Ag": 2,
    "CD": 2,
    "PD": 2,
    "PT": 2,
    "HG": 2,
    "AL": 3,
    "IN": 3,
    "TL": 1,
    "SN": 2,
    "Sn": 2,
    "PB": 2,
    "PR": 3,
    "ND": 3,
    "SM": 3,
    "Sm": 2,
    "EU": 2,
    "EU3": 3,
    "GD3": 3,
    "TB": 3,
    "Dy": 3,
    "Er": 3,
    "Tm": 3,
    "YB2": 2,
    # Actinides
    "Th": 4,
    "U4+": 4,
    "Pu": 3,
    # Halogens
    "F": -1,
    "CL": -1,
    "Cl-": -1,
    "BR": -1,
    "IOD": -1,
    # Other common ions
    "H3O+": 1,
    "NH4": 1,
    "HZ+": 1,
    "HE+": 1,
    # Other metals
    "Zr": 4,
    "Hf": 4,
}


class ProteinComponent(ExplicitMoleculeComponent):
    """
    :class:`Component` representing the contents of a PDB file, such as a protein.

    In comparison to a SmallMoleculeComponent, this representation additionally
    contains information relating to the residue and chain information.  Technically,
    this is done by having the ``MonomerInfo`` attributes present on each atom
    of the input RDKit molecule, which is done when reading from either PDB or
    ``.mae`` file inputs.

    Parameters
    ----------
    rdkit : rdkit.Mol
       rdkit representation of the protein
    name : str, optional
       of the protein, by default ""

    Note
    ----
    This class is a read-only representation of a protein, if you want to
    edit the molecule do this in an appropriate toolkit **before** creating
    an instance from this class.
    """

    def __init__(self, rdkit: RDKitMol, name=""):
        if not all(a.GetMonomerInfo() is not None for a in rdkit.GetAtoms()):
            raise TypeError(
                "Not all atoms in input have MonomerInfo defined.  "
                "Consider loading via rdkit.Chem.MolFromPDBFile or similar."
            )
        super().__init__(rdkit=rdkit, name=name)

    # FROM
    @classmethod
    def from_pdb_file(cls, pdb_file: str, name: str = ""):
        """
        Create ``ProteinComponent`` from PDB-formatted file.

        Parameters
        ----------
        pdb_file : str
            path to the pdb file.
        name : str, optional
            name of the input protein, by default ""

        Returns
        -------
        ProteinComponent
            the deserialized molecule
        """
        openmm_PDBFile = PDBFile(pdb_file)
        return cls._from_openmmPDBFile(openmm_PDBFile=openmm_PDBFile, name=name)

    @classmethod
    def from_pdbx_file(cls, pdbx_file: str, name=""):
        """
        Create ``ProteinComponent`` from PDBX-formatted file.

        Parameters
        ----------
        pdbxfile : str
            path to the pdb file.
        name : str, optional
            name of the input protein, by default ""

        Returns
        -------
        ProteinComponent
            the deserialized molecule
        """
        openmm_PDBxFile = PDBxFile(pdbx_file)
        return cls._from_openmmPDBFile(openmm_PDBFile=openmm_PDBxFile, name=name)

    @classmethod
    def _from_openmmPDBFile(cls, openmm_PDBFile: PDBFile | PDBxFile, name: str = ""):
        """Converts to our internal representation (rdkit Mol)

        Parameters
        ----------
        openmm_PDBFile : PDBFile or PDBxFile
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

        # Add Atoms
        for atom in mol_topology.atoms():
            a = Atom(atom.element.atomic_number)

            atom_monomerInfo = Chem.AtomPDBResidueInfo()
            atom_monomerInfo.SetChainId(atom.residue.chain.id)
            atom_monomerInfo.SetSerialNumber(int(atom.id))
            atom_monomerInfo.SetSegmentNumber(int(atom.residue.chain.index))
            atom_monomerInfo.SetInsertionCode(atom.residue.insertionCode)
            atom_monomerInfo.SetName(atom.name)
            atom_monomerInfo.SetResidueName(atom.residue.name)
            atom_monomerInfo.SetResidueNumber(int(atom.residue.id))
            atom_monomerInfo.SetIsHeteroAtom(False)  # TODO: Do hetatoms

            a.SetMonomerInfo(atom_monomerInfo)

            # additonally possible:
            # atom_monomerInfo.SetSecondaryStructure
            # atom_monomerInfo.SetMonomerType
            # atom_monomerInfo.SetAltLoc

            editable_rdmol.AddAtom(a)

        # Add Bonds
        for bond in mol_topology.bonds():
            bond_order = _BONDORDERS_OPENMM_TO_RDKIT[bond.order]
            editable_rdmol.AddBond(
                beginAtomIdx=bond.atom1.index,
                endAtomIdx=bond.atom2.index,
                order=bond_order,
            )

        # Set Positions
        rd_mol = editable_rdmol.GetMol()
        positions = np.array(openmm_PDBFile.positions.value_in_unit(omm_unit.angstrom), ndmin=3)

        for frame_id, frame in enumerate(positions):
            conf = Conformer(frame_id)
            for atom_id, atom_pos in enumerate(frame):
                conf.SetAtomPosition(atom_id, atom_pos)
            rd_mol.AddConformer(conf)

        # Add Additionals
        # Formal Charge
        netcharge = 0
        for a in rd_mol.GetAtoms():
            atom_name = a.GetMonomerInfo().GetName().strip()
            atomic_num = a.GetAtomicNum()

            connectivity = sum(_BONDORDER_TO_ORDER[bond.GetBondType()] for bond in a.GetBonds())
            default_valence = periodicTable.GetDefaultValence(atomic_num)

            if connectivity == 0:  # ions
                ion_key = atom_name.strip().upper()
                if ion_key in ions_dict:
                    fc = ions_dict[ion_key]
                else:
                    resn = a.GetMonomerInfo().GetResidueName()
                    resind = int(a.GetMonomerInfo().GetResidueNumber())
                    raise ValueError(
                        f"Unknown ion: {atom_name} in residue {resn} at index {resind}. "
                        f"Check if it's in the ions_dict dictionary."
                    )
            elif default_valence > connectivity:
                fc = -(default_valence - connectivity)  # negative charge
            elif default_valence < connectivity:
                fc = +(connectivity - default_valence)  # positive charge
            else:
                fc = 0  # neutral

            a.SetFormalCharge(fc)
            a.UpdatePropertyCache(strict=True)

            netcharge += fc

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

            a = Atom(atomic_num)
            mi = Chem.AtomPDBResidueInfo()

            mi.SetChainId(atom[1])
            mi.SetSerialNumber(atom[2])
            mi.SetSegmentNumber(atom[3])
            mi.SetInsertionCode(atom[4])
            mi.SetName(atom[5])
            mi.SetResidueName(atom[6])
            mi.SetResidueNumber(int(atom[7]))
            mi.SetIsHeteroAtom(atom[8] == "Y")
            a.SetFormalCharge(atom[9])

            a.SetMonomerInfo(mi)

            editable_rdmol.AddAtom(a)

        # Add Bonds
        for bond in ser_dict["bonds"]:
            atomBeginIdx = int(bond[0])
            atomEndIdx = int(bond[1])
            bondType = _BONDORDER_STR_TO_RDKIT[bond[2]]
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
            # Can't set these on an editable mol, go round a second time
            _, _, _, arom = ser_dict["bonds"][bond_id]
            bond.SetIsAromatic(arom == "Y")

        if "name" in ser_dict:
            name = ser_dict["name"]

        return cls(rdkit=rd_mol, name=name)

    def to_openmm_topology(self) -> app.Topology:
        """Convert to an openmm Topology object

        Returns
        -------
        openmm.app.Topology
            resulting topology obj.
        """

        def reskey(m):
            """key for defining when a residue has changed from previous

            this matches criteria used in openmm (pdbstructure), except altloc
            ignored
            """
            return (
                m.GetChainId(),
                m.GetResidueName(),
                m.GetResidueNumber(),
                m.GetInsertionCode(),
            )

        def chainkey(m):
            """key for chains

            uses (chain.id, chain.index)

            using .index catches where TER records in openmm have caused a new
            chain to be used (with identical .id)
            """
            return (
                m.GetChainId(),
                m.GetSegmentNumber(),
            )

        current_chainid = None
        c = None  # current chain
        current_resid = None
        r = None  # current residue

        atom_lookup = {}  # maps rdkit indices to openmm Atoms

        top = app.Topology()
        for atom in self._rdkit.GetAtoms():
            mi = atom.GetMonomerInfo()
            if (new_chainid := chainkey(mi)) != current_chainid:
                chainid, _ = new_chainid
                c = top.addChain(chainid)
                current_chainid = new_chainid

            if (new_resid := reskey(mi)) != current_resid:
                _, resname, resnum, icode = new_resid
                r = top.addResidue(name=resname, chain=c, id=str(resnum), insertionCode=icode)
                current_resid = new_resid

            a = top.addAtom(
                name=mi.GetName(),
                element=app.Element.getByAtomicNumber(atom.GetAtomicNum()),
                residue=r,
                id=str(mi.GetSerialNumber()),
            )

            atom_lookup[atom.GetIdx()] = a

        for bond in self._rdkit.GetBonds():
            a1 = atom_lookup[bond.GetBeginAtomIdx()]
            a2 = atom_lookup[bond.GetEndAtomIdx()]
            rdkit_bond_type = bond.GetBondType()
            bond_order = _BONDORDERS_RDKIT_TO_OPENMM.get(rdkit_bond_type, None)
            bond_type = _BONDTYPES_RDKIT_TO_OPENMM.get(rdkit_bond_type, None)
            top.addBond(a1, a2, order=bond_order, type=bond_type)

        return top

    def to_openmm_positions(self) -> omm_unit.Quantity:
        """
        serialize the positions to openmm.unit.Quantity

        Note
        ----
        Currently only one frame/model is given

        Returns
        -------
        omm_unit.Quantity
            Quantity containing protein atom positions
        """
        np_pos = deserialize_numpy(self.to_dict()["conformers"][0])
        openmm_pos = list(map(lambda x: np.array(x), np_pos)) * omm_unit.angstrom

        return openmm_pos

    def to_pdb_file(self, out_path: str | bytes | PathLike[str] | PathLike[bytes] | io.TextIOBase) -> str:
        """
        serialize protein to pdb file.

        Parameters
        ----------
        out_path :  Union[str, bytes, PathLike[str], PathLike[bytes], io.TextIOBase]
            provide path or any string based stream (e.g. FileIO ) to the resulting file

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
        if not isinstance(out_path, io.TextIOBase):
            # allows pathlike/str; we close on completion
            out_file = open(out_path, mode="w")  # type: ignore
            must_close = True
        else:
            out_file = out_path  # type: ignore
            must_close = False

        try:
            out_path = out_file.name
        except AttributeError:
            out_path = "<unknown>"

        PDBFile.writeFile(topology=openmm_top, positions=openmm_pos, file=out_file)

        if must_close:
            # we only close the file if we had to open it
            out_file.close()

        return out_path

    def to_pdbx_file(self, out_path: str | bytes | PathLike[str] | PathLike[bytes] | io.TextIOBase) -> str:
        """
        serialize protein to pdbx file.

        Parameters
        ----------
        out_path : Union[str, bytes, PathLike[str], PathLike[bytes], io.TextIOBase]
            provide path or FileIO to the resulting file

        Returns
        -------
        str
            string path to the resulting pdbx.
        """
        # get top:
        top = self.to_openmm_topology()

        # get pos:
        np_pos = deserialize_numpy(self.to_dict()["conformers"][0])
        openmm_pos = list(map(lambda x: np.array(x), np_pos)) * omm_unit.angstrom

        # write file
        if not isinstance(out_path, io.TextIOBase):
            # allows pathlike/str; we close on completion
            out_file = open(out_path, mode="w")  # type: ignore
            must_close = True
        else:
            out_file = out_path  # type: ignore
            must_close = False

        try:
            out_path = out_file.name
        except AttributeError:
            out_path = "<unknown>"

        PDBxFile.writeFile(topology=top, positions=openmm_pos, file=out_file)

        if must_close:
            # we only close the file if we had to open it
            out_file.close()

        return out_path

    def _to_dict(self) -> dict:
        """Serialize to dict representation"""

        atoms = []
        for atom in self._rdkit.GetAtoms():
            mi = atom.GetMonomerInfo()

            # TODO: Stereo?
            atoms.append(
                (
                    atom.GetAtomicNum(),
                    mi.GetChainId(),
                    mi.GetSerialNumber(),
                    mi.GetSegmentNumber(),
                    mi.GetInsertionCode(),
                    mi.GetName(),
                    mi.GetResidueName(),
                    mi.GetResidueNumber(),
                    "Y" if mi.GetIsHeteroAtom() else "N",
                    atom.GetFormalCharge(),
                )
            )

        bonds = [
            (
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                _BONDORDER_RDKIT_TO_STR[bond.GetBondType()],
                "Y" if bond.GetIsAromatic() else "N",
                # bond.GetStereo() or "",  do we need this? i.e. are openff ffs going to use cis/trans SMARTS?
            )
            for bond in self._rdkit.GetBonds()
        ]

        conformers = [
            serialize_numpy(conf.GetPositions()) for conf in self._rdkit.GetConformers()  # .m_as(unit.angstrom)
        ]

        # Result
        d = {
            "atoms": atoms,
            "bonds": bonds,
            "name": self.name,
            "conformers": conformers,
        }

        return d
