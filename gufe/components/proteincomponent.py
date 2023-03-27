# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import io
import numpy as np
from os import PathLike
from typing import Union, Optional
import pdbinf

from rdkit import Chem

from ..custom_typing import RDKitMol
from .explicitmoleculecomponent import ExplicitMoleculeComponent, _mol_from_dict


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
    """Wrapper around a Protein representation.

    In comparison to a SmallMoleculeComponent, this representation additionally contains information
    relating to the residue and chain information.  This is achievable by having the MonomerInfo attributes
    present on each atom of the input RDKit molecule, which is done when reading from either PDB or `.mae`
    file inputs.

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
        """Create ``ProteinComponent`` from PDB-formatted file.

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
        return cls(rdkit=pdbinf.load_pdb_file(pdb_file,
                                              templates=[pdbinf.STANDARD_AA_DOC, pdbinf.DNA_DOC, pdbinf.RNA_DOC],
                                              ),
                   name=name)

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
        return cls(rdkit=pdbinf.load_pdbx_file(pdbx_file,
                                               templates=[pdbinf.STANDARD_AA_DOC, pdbinf.DNA_DOC, pdbinf.RNA_DOC],
                                               ),
                   name=name)

    @classmethod
    def _from_dict(cls, ser_dict: dict):
        """Deserialize from dict representation"""
        nm = ser_dict.pop('name', '')
        monomerinfo = ser_dict.pop('monomerinfo')

        m = _mol_from_dict(ser_dict)

        for at, mi in zip(m.GetAtoms(), monomerinfo):
            info = Chem.AtomPDBResidueInfo()
            info.SetSerialNumber(mi[0])
            info.SetAltLoc(mi[1])
            info.SetResidueName(mi[2])
            info.SetResidueNumber(mi[3])
            info.SetChainId(mi[4])
            info.SetInsertionCode(mi[5])
            info.SetOccupancy(mi[6])
            info.SetTempFactor(mi[7])
            info.SetIsHeteroAtom(mi[8])
            info.SetSecondaryStructure(mi[9])
            info.SetSegmentNumber(mi[10])

            at.SetMonomerInfo(info)

        return cls(rdkit=m, name=nm)

    def to_openmm_topology(self):
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
                r = top.addResidue(name=resname.strip(),
                                   chain=c,
                                   id=str(resnum),
                                   insertionCode=icode)
                current_resid = new_resid

            a = top.addAtom(
                name=mi.GetName().strip(),  # PDBFile seems to strip whitespace, rdkit doesn't
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

    def to_openmm_positions(self):
        """Convert the positions to openmm.unit.Quantity

        ! only one frame at the moment!

        Returns
        -------
        omm_unit.Quantity
            Quantity containing protein atom positions
        """
        from openmm import unit as omm_unit

        np_pos = self._rdkit.GetConformer().GetPositions()

        return np_pos * omm_unit.angstrom

    def to_pdb_file(self, out_path: Union[str, bytes, PathLike[str], PathLike[bytes], io.TextIOBase]) -> str:
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
        from openmm.app import PDBFile

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
        from openmm.app import PDBxFile

        # get top:
        top = self.to_openmm_topology()
        pos = self.to_openmm_positions()

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

        PDBxFile.writeFile(topology=top, positions=pos, file=out_file)

        if must_close:
            # we only close the file if we had to open it
            out_file.close()

        return out_path

    def _to_dict(self) -> dict:
        """Serialize to dict representation"""
        # ExplicitMolecule does everything except monomer info,
        # so use that then augment with MonomerInfo
        d = super()._to_dict()

        monomerinfo = []
        for atom in self._rdkit.GetAtoms():
            mi: Chem.AtomPDBResidueInfo = atom.GetMonomerInfo()
            monomerinfo.append((
                mi.GetSerialNumber(), mi.GetAltLoc(), mi.GetResidueName(), mi.GetResidueNumber(), mi.GetChainId(),
                mi.GetInsertionCode(), mi.GetOccupancy(), mi.GetTempFactor(), mi.GetIsHeteroAtom(),
                mi.GetSecondaryStructure(), mi.GetSegmentNumber()
            ))
            # monomerinfo structs cannot carry props

        d['monomerinfo'] = monomerinfo

        return d
