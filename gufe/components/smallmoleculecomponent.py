# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import logging
# openff complains about oechem being missing, shhh
logger = logging.getLogger('openff.toolkit')
logger.setLevel(logging.ERROR)
from typing import Any

from rdkit import Chem

from .explicitmoleculecomponent import ExplicitMoleculeComponent
from ..molhashing import deserialize_numpy, serialize_numpy


_INT_TO_ATOMCHIRAL = {
    0: Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
    1: Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
    2: Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW,
    3: Chem.rdchem.ChiralType.CHI_OTHER,
}
# support for non-tetrahedral stereo requires rdkit 2022.09.1+
if hasattr(Chem.rdchem.ChiralType, 'CHI_TETRAHEDRAL'):
    _INT_TO_ATOMCHIRAL.update({
        4: Chem.rdchem.ChiralType.CHI_TETRAHEDRAL,
        5: Chem.rdchem.ChiralType.CHI_ALLENE,
        6: Chem.rdchem.ChiralType.CHI_SQUAREPLANAR,
        7: Chem.rdchem.ChiralType.CHI_TRIGONALBIPYRAMIDAL,
        8: Chem.rdchem.ChiralType.CHI_OCTAHEDRAL,
    })
_ATOMCHIRAL_TO_INT = {v: k for k, v in _INT_TO_ATOMCHIRAL.items()}


_INT_TO_BONDTYPE = {
    0: Chem.rdchem.BondType.UNSPECIFIED,
    1: Chem.rdchem.BondType.SINGLE,
    2: Chem.rdchem.BondType.DOUBLE,
    3: Chem.rdchem.BondType.TRIPLE,
    4: Chem.rdchem.BondType.QUADRUPLE,
    5: Chem.rdchem.BondType.QUINTUPLE,
    6: Chem.rdchem.BondType.HEXTUPLE,
    7: Chem.rdchem.BondType.ONEANDAHALF,
    8: Chem.rdchem.BondType.TWOANDAHALF,
    9: Chem.rdchem.BondType.THREEANDAHALF,
    10: Chem.rdchem.BondType.FOURANDAHALF,
    11: Chem.rdchem.BondType.FIVEANDAHALF,
    12: Chem.rdchem.BondType.AROMATIC,
    13: Chem.rdchem.BondType.IONIC,
    14: Chem.rdchem.BondType.HYDROGEN,
    15: Chem.rdchem.BondType.THREECENTER,
    16: Chem.rdchem.BondType.DATIVEONE,
    17: Chem.rdchem.BondType.DATIVE,
    18: Chem.rdchem.BondType.DATIVEL,
    19: Chem.rdchem.BondType.DATIVER,
    20: Chem.rdchem.BondType.OTHER,
    21: Chem.rdchem.BondType.ZERO}
_BONDTYPE_TO_INT = {v: k for k, v in _INT_TO_BONDTYPE.items()}
_INT_TO_BONDSTEREO = {
    0: Chem.rdchem.BondStereo.STEREONONE,
    1: Chem.rdchem.BondStereo.STEREOANY,
    2: Chem.rdchem.BondStereo.STEREOZ,
    3: Chem.rdchem.BondStereo.STEREOE,
    4: Chem.rdchem.BondStereo.STEREOCIS,
    5: Chem.rdchem.BondStereo.STEREOTRANS}
_BONDSTEREO_TO_INT = {v: k for k, v in _INT_TO_BONDSTEREO.items()}


def _setprops(obj, d: dict) -> None:
    # add props onto rdkit "obj" (atom/bond/mol/conformer)
    # props are guaranteed one of Bool, Int, Float or String type
    for k, v in d.items():
        if isinstance(v, bool):
            obj.SetBoolProp(k, v)
        elif isinstance(v, int):
            obj.SetIntProp(k, v)
        elif isinstance(v, float):
            obj.SetDoubleProp(k, v)
        else:  # isinstance(v, str):
            obj.SetProp(k, v)


class SmallMoleculeComponent(ExplicitMoleculeComponent):
    """
    :class:`Component` representing a small molecule, used for ligands and cofactors.

    This class supports reading/writing to the ``.sdf`` format, which is suited to
    smaller molecules. Ligands in a free energy calculation are represented with this class,
    along with cofactors.

    The name can be explicitly set by the ``name`` keyword argument on create,
    or implicitly set based on the tags in the input molecular representation
    (if supported, as with RDKit). If not explicitly set on creation, the molecule
    will first look for an OpenFE-specific tag ``ofe-name``, and if that doesn't exist,
    for a commonly-used naming tag (e.g., the ``_Name`` property for RDKit
    molecules). If no name is found, the empty string is used.

    Parameters
    ----------
    rdkit : |rdkit.mol|
        rdkit representation of the molecule
    name : str, optional
        A human readable tag for this molecule.  This name will be used in the hash.

    Note
    ----
    This class is a read-only representation of a molecule, if you want to
    edit the molecule do this in an appropriate toolkit **before** creating
    an instance from this class.
    """

    def to_sdf(self) -> str:
        """Create a string based on SDF.

        See Also
        --------
        :meth:`.from_sdf_string` : create an object from the output of this
        """
        # https://sourceforge.net/p/rdkit/mailman/message/27518272/
        mol = self.to_rdkit()
        sdf = [Chem.MolToMolBlock(mol)]
        for prop in mol.GetPropNames():
            val = mol.GetProp(prop)
            sdf.append('>  <%s>\n%s\n' % (prop, val))
        sdf.append('$$$$\n')
        return "\n".join(sdf)

    @classmethod
    def from_sdf_string(cls, sdf_str: str):
        """Create `:class:SmallMoleculeComponent` from SDF-formatted string.

        Parameters
        ----------
        sdf_str : str
            input string in SDF format

        Returns
        -------
        :class:`.SmallMoleculeComponent` :
            the deserialized molecule
        """
        supp = Chem.SDMolSupplier()
        supp.SetData(sdf_str, removeHs=False)
        return cls._from_sdf_supplier(supp)

    @classmethod
    def from_sdf_file(cls, filename: str):
        """Create ``SmallMoleculeComponent`` from SDF file.

        Parameters
        ----------
        filename : str
            name of SDF file

        Returns
        -------
        :class:`.SmallMoleculeComponent` :
            the deserialized molecule
        """
        # technically, we allow file-like objects
        supp = Chem.SDMolSupplier(str(filename), removeHs=False)
        return cls._from_sdf_supplier(supp)

    @classmethod
    def _from_sdf_supplier(cls, supp):
        """
        Internal mechanism used by both from_sdf_string and from_sdf_file.
        """
        mol = next(supp)
        if mol is None:
            raise ValueError("Unable to load SmallMoleculeComponent")

        # ensure that there's only one molecule in the file
        try:
            _ = next(supp)
        except StopIteration:
            pass
        else:
            # TODO: less generic exception type here
            raise RuntimeError(f"SDF contains more than 1 molecule")

        return cls(rdkit=mol)  # name is obtained automatically

    def to_openff(self):
        """OpenFF Toolkit Molecule representation of this molecule

        Note
        ----
        This is a copy of this object, and modifying the OpenFF copy does not
        alter the original object.
        """
        from openff.toolkit.topology import Molecule as OFFMolecule

        m = OFFMolecule(self._rdkit, allow_undefined_stereo=True)
        m.name = self.name

        return m

    @classmethod
    def from_openff(cls, openff, name: str = ""):
        """Construct from an OpenFF toolkit Molecule"""
        return cls(openff.to_rdkit(), name=name)

    def _to_dict(self) -> dict:
        """Serialize to dict representation"""
        # in a perfect world we'd use ToBinary()
        # but this format slowly evolves, so the future hash of a SMC could change if rdkit were updated
        # this is based on that method, with some irrelevant fields cut out

        output: dict[str, Any] = {}

        atoms = []
        for atom in self._rdkit.GetAtoms():
            atoms.append((
                atom.GetAtomicNum(), atom.GetIsotope(), atom.GetFormalCharge(), atom.GetIsAromatic(),
                _ATOMCHIRAL_TO_INT[atom.GetChiralTag()], atom.GetAtomMapNum(),
                atom.GetPropsAsDict(includePrivate=False),
            ))
        output['atoms'] = atoms

        bonds = []
        for bond in self._rdkit.GetBonds():
            bonds.append((
                bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), _BONDTYPE_TO_INT[bond.GetBondType()],
                _BONDSTEREO_TO_INT[bond.GetStereo()],
                bond.GetPropsAsDict(includePrivate=False)
            ))
        output['bonds'] = bonds

        conf = self._rdkit.GetConformer()
        output['conformer'] = (serialize_numpy(conf.GetPositions()), conf.GetPropsAsDict(includePrivate=False))

        output['molprops'] = self._rdkit.GetPropsAsDict(includePrivate=False)

        return output

    @classmethod
    def _from_dict(cls, d: dict):
        """Deserialize from dict representation"""
        m = Chem.Mol()
        em = Chem.EditableMol(m)

        for atom in d['atoms']:
            a = Chem.Atom(atom[0])
            a.SetIsotope(atom[1])
            a.SetFormalCharge(atom[2])
            a.SetIsAromatic(atom[3])
            a.SetChiralTag(_INT_TO_ATOMCHIRAL[atom[4]])
            a.SetAtomMapNum(atom[5])
            _setprops(a, atom[6])
            em.AddAtom(a)

        for bond in d['bonds']:
            em.AddBond(bond[0], bond[1], _INT_TO_BONDTYPE[bond[2]])
            # other fields are applied onto the ROMol

        m = em.GetMol()

        for bond, b in zip(d['bonds'], m.GetBonds()):
            b.SetStereo(_INT_TO_BONDSTEREO[bond[3]])
            _setprops(b, bond[4])

        pos = deserialize_numpy(d['conformer'][0])
        c = Chem.Conformer(m.GetNumAtoms())
        for i, p in enumerate(pos):
            c.SetAtomPosition(i, p)
        _setprops(c, d['conformer'][1])
        m.AddConformer(c)

        _setprops(m, d['molprops'])

        m.UpdatePropertyCache()

        return cls(rdkit=m)

    def copy_with_replacements(self, **replacements):
        # this implementation first makes a copy with the name replaced
        # only, then does any other replacements that are necessary
        if 'name' in replacements:
            name = replacements.pop('name')
            dct = self._to_dict()
            dct['molprops']['ofe-name'] = name
            obj = self._from_dict(dct)
        else:
            obj = self

        return super(SmallMoleculeComponent, obj).copy_with_replacements(**replacements)
