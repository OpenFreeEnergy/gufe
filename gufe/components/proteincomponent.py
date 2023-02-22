# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import io
import numpy as np
from os import PathLike
from typing import Union

import pdbinf

from rdkit import Chem
from rdkit.Chem.rdchem import Mol, Atom, Conformer, EditableMol, BondType

from ..custom_typing import RDKitMol
from .explicitmoleculecomponent import ExplicitMoleculeComponent
from ..vendor.pdb_file.pdbfile import PDBFile
from ..vendor.pdb_file.pdbxfile import PDBxFile


from ..molhashing import deserialize_numpy, serialize_numpy


_BONDORDERS_OPENMM_TO_RDKIT = {
    1: BondType.SINGLE,
    2: BondType.DOUBLE,
    3: BondType.TRIPLE,
    app.Single: BondType.SINGLE,
    app.Double: BondType.DOUBLE,
    app.Triple: BondType.TRIPLE,
    app.Aromatic: BondType.AROMATIC,
    None: BondType.UNSPECIFIED,
}
_BONDORDERS_RDKIT_TO_OPENMM = {
    v: k for k, v in _BONDORDERS_OPENMM_TO_RDKIT.items()
}

# builtin dict of strings to enum members, boy I hope this is stable
_BONDORDER_STR_TO_RDKIT = Chem.BondType.names
_BONDORDER_RDKIT_TO_STR = {v: k for k, v in _BONDORDER_STR_TO_RDKIT.items()}

_CHIRALITY_RDKIT_TO_STR = {
    Chem.CHI_TETRAHEDRAL_CW: 'CW',
    Chem.CHI_TETRAHEDRAL_CCW: 'CCW',
    Chem.CHI_UNSPECIFIED: 'U',
}
_CHIRALITY_STR_TO_RDKIT = {
    v: k for k, v in _CHIRALITY_RDKIT_TO_STR.items()
}


class ProteinComponent(ExplicitMoleculeComponent):
    """
    ``Component`` representing the contents of a PDB file, such as a protein.

    In comparison to a SmallMoleculeComponent, this representation additionally
    contains information relating to the residue and chain information.  This
    is achievable by having the ``MonomerInfo`` attributes present on each atom
    of the input RDKit molecule, which is done when reading from either PDB or
    ``.mae`` file inputs.

    Note
    ----
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
    def __init__(self, rdkit: RDKitMol, name=""):
        if not all(a.GetMonomerInfo() is not None for a in rdkit.GetAtoms()):
            raise TypeError("Not all atoms in input have MonomerInfo defined.  "
                            "Consider loading via rdkit.Chem.MolFromPDBFile or similar.")
        super().__init__(rdkit=rdkit, name=name)

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
        return cls(
            rdkit=pdbinf.load_pdb_file(pdb_file,
                                 templates=[pdbinf.STANDARD_AA_DOC]),
            name=name,
        )

    @classmethod
    def from_pdbx_file(cls, pdbx_file: str, name=""):
        """
        Create ``ProteinComponent`` from PDBX-formatted file.

        Parameters
        ----------
        pdbx_file : str
            path to the pdb file.
        name : str, optional
            name of the input protein, by default ""

        Returns
        -------
        ProteinComponent
            the deserialized molecule
        """
        return cls(
            rdkit=pdbinf.load_pdbx_file(pdbx_file,
                                  templates=[pdbinf.STANDARD_AA_DOC]),
            name=name,
        )

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
            mi.SetIsHeteroAtom(atom[8] == 'Y')
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
            bond.SetIsAromatic(arom == 'Y')

        if "name" in ser_dict:
            name = ser_dict["name"]

        return cls(rdkit=rd_mol, name=name)

    def to_openmm_topology(self) -> "app.Topology":
        """Convert to an openmm Topology object

        Returns
        -------
        openmm.app.Topology
            resulting topology obj.
        """
        from openmm import app

        def reskey(m):
            """key for defining when a residue has changed from previous

            this matches criteria used in openmm (pdbstructure), except altloc
            ignored
            """
            return (
                m.GetChainId(),
                m.GetResidueName(),
                m.GetResidueNumber(),
                m.GetInsertionCode()
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
                r = top.addResidue(name=resname,
                                   chain=c,
                                   id=str(resnum),
                                   insertionCode=icode)
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
            top.addBond(a1, a2,
                        order=_BONDORDERS_RDKIT_TO_OPENMM.get(
                            bond.GetBondType(), None))

        return top

    def to_openmm_positions(self) -> "openmm.app.unit.Quantity":
        """serialize the positions to openmm.unit.Quantity

        ! only one frame at the moment!

        Returns
        -------
        omm_unit.Quantity
            Quantity containing protein atom positions
        """
        from openmm import unit as omm_unit

        np_pos = deserialize_numpy(self.to_dict()["conformers"][0])
        openmm_pos = (
            list(map(lambda x: np.array(x), np_pos)) * omm_unit.angstrom
        )

        return openmm_pos

    def to_pdb_file(self, out_path: Union[str, bytes, PathLike[str], PathLike[bytes], io.TextIOBase]) -> str:
        """Write protein to pdb file.

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
            out_file = open(out_path, mode='w')  # type: ignore
            must_close = True
        else:
            out_file = out_path  # type: ignore
            must_close = False
            
        try:
            out_path = out_file.name
        except AttributeError:
            out_path = "<unknown>"

        PDBFile.writeFile(
            topology=openmm_top, positions=openmm_pos, file=out_file
        )

        if must_close:
            # we only close the file if we had to open it
            out_file.close()

        return out_path

    def to_pdbx_file(
        self, out_path: Union[str, bytes, PathLike[str], PathLike[bytes], io.TextIOBase]
    ) -> str:
        """Write protein to pdbx file.

        Parameters
        ----------
        out_path : Union[str, bytes, PathLike[str], PathLike[bytes], io.TextIOBase]
            provide path or FileIO to the resulting file

        Returns
        -------
        str
            string path to the resulting pdbx.
        """
        from openmm import unit as omm_unit

        # get top:
        top = self.to_openmm_topology()

        # get pos:
        np_pos = deserialize_numpy(self.to_dict()["conformers"][0])
        openmm_pos = (
            list(map(lambda x: np.array(x), np_pos)) * omm_unit.angstrom
        )

        # write file
        if not isinstance(out_path, io.TextIOBase):
            # allows pathlike/str; we close on completion
            out_file = open(out_path, mode='w')  # type: ignore
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
                    'Y' if mi.GetIsHeteroAtom() else 'N',
                    atom.GetFormalCharge(),
                )
            )

        bonds = [
            (
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                _BONDORDER_RDKIT_TO_STR[bond.GetBondType()],
                'Y' if bond.GetIsAromatic() else 'N',
                # bond.GetStereo() or "",  do we need this? i.e. are openff ffs going to use cis/trans SMARTS?
            )
            for bond in self._rdkit.GetBonds()
        ]

        conformers = [
            serialize_numpy(conf.GetPositions())  # .m_as(unit.angstrom)
            for conf in self._rdkit.GetConformers()
        ]

        # Result
        d = {
            "atoms": atoms,
            "bonds": bonds,
            "name": self.name,
            "conformers": conformers,
        }

        return d
