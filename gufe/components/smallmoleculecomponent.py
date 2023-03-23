# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import logging
# openff complains about oechem being missing, shhh
logger = logging.getLogger('openff.toolkit')
logger.setLevel(logging.ERROR)

from rdkit import Chem

from .explicitmoleculecomponent import ExplicitMoleculeComponent, _mol_from_dict

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
    """A molecule wrapper suitable for small molecules

    .. note::
       This class is a read-only representation of a molecule, if you want to
       edit the molecule do this in an appropriate toolkit **before** creating
       an instance from this class.

    This class supports reading/writing to the `.sdf` format, which is suited to smaller
    molecules.

    The name can be explicitly set by the ``name`` attribute, or implicitly set
    based on the tags in the input molecular representation (if supported, as
    with RDKit). If not explicitly set on creation, the molecule will first
    look for an OpenFE-specific tag ``ofe-name``, and if that doesn't exist,
    for a commonly-used naming tag (e.g., the ``_Name`` property for RDKit
    molecules). If no name is found, the empty string is used.

    Parameters
    ----------
    rdkit : |rdkit.mol|
        rdkit representation of the molecule
    name : str, optional
        A human readable tag for this molecule.  This name will be used in the hash.
    """
    @classmethod
    def _from_dict(cls, dct: dict):
        nm = dct.pop('name', '')
        m = _mol_from_dict(dct)

        return cls(rdkit=m, name=nm)

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
        """Create ``SmallMoleculeComponent`` from SDF-formatted string.

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
        """OpenFF Toolkit representation of this molecule"""
        from openff.toolkit.topology import Molecule as OFFMolecule

        m = OFFMolecule(self._rdkit, allow_undefined_stereo=True)
        m.name = self.name

        return m

    @classmethod
    def from_openff(cls, openff, name: str = ""):
        """Construct from an OpenFF toolkit Molecule"""
        return cls(openff.to_rdkit(), name=name)
