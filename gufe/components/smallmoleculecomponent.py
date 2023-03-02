# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import logging
# openff complains about oechem being missing, shhh
logger = logging.getLogger('openff.toolkit')
logger.setLevel(logging.ERROR)
from openff.toolkit.topology import Molecule as OFFMolecule
from openff.units import unit

from rdkit import Chem

from .explicitmoleculecomponent import ExplicitMoleculeComponent
from ..custom_typing import OEMol
from ..molhashing import deserialize_numpy, serialize_numpy



class SmallMoleculeComponent(ExplicitMoleculeComponent):
    """A molecule wrapper suitable for small molecules

    .. note::
       This class is a read-only representation of a molecule, if you want to
       edit the molecule do this in an appropriate toolkit **before** creating
       an instance from this class.

    This wrapper uses SMILES as the primary hash, so is best suited to smaller
    molecules.  It also supports reading/writing to .sdf format, which again
    is suited to smaller molecules.
    A small molecule can have a name associated with it, which is needed to
    distinguish two molecules with the same SMILES representation, or is
    simply useful to help identify ligand molecules later.
    The name can be explicitly set by the ``name`` attribute, or implicitly set
    based on the tags in the input molecular representation (if supported, as
    with RDKit). If not explicitly set on creation, the molecule will first
    look for an OpenFE-specific tag ``ofe-name``, and if that doesn't exist,
    for a commonly-used naming tag (e.g., the ``_Name`` property for RDKit
    molecules). If no name is found, the empty string is used.

    Parameters
    ----------
    rdkit : rdkit.Mol
        rdkit representation of the molecule
    name : str, optional
        if multiple Molecules with identical SMILES but differing positions
        are used, a name must be given to differentiate these.  This name
        will be used in the hash.
    """

    def to_sdf(self) -> str:
        """Create a string based on SDF.

        This is the primary serialization mechanism for this class.

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

        This is the primary deserialization mechanism for this class.

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

    def to_openeye(self) -> OEMol:  # type: ignore
        # typing: see https://github.com/OpenFreeEnergy/gufe/pull/65#issuecomment-1259682099
        """OEChem representation of this molecule"""
        return self.to_openff().to_openeye()

    @classmethod
    def from_openeye(cls, oemol: OEMol, name: str = ""):
        raise NotImplementedError

    def to_openff(self):
        """OpenFF Toolkit representation of this molecule"""
        m = OFFMolecule(self._rdkit, allow_undefined_stereo=True)
        m.name = self.name

        return m

    @classmethod
    def from_openff(cls, openff: OFFMolecule, name: str = ""):
        """Construct from an OpenFF toolkit Molecule"""
        return cls(openff.to_rdkit(), name=name)

    def _to_dict(self) -> dict:
        """Serialize to dict representation"""
        # this changes a global property
        before = Chem.GetDefaultPickleProperties()

        Chem.SetDefaultPickleProperties(
            Chem.PropertyPickleOptions.AtomProps |
            Chem.PropertyPickleOptions.MolProps |
            Chem.PropertyPickleOptions.BondProps |
            Chem.PropertyPickleOptions.CoordsAsDouble
        )

        blob = self._rdkit.ToBinary()

        Chem.SetDefaultPickleProperties(before)

        return {
            'rdkit_blob': blob,
            'name': self._name,
        }

    @classmethod
    def _from_dict(cls, d: dict):
        """Deserialize from dict representation"""
        m = Chem.Mol(d['rdkit_blob'])

        return cls(rdkit=m, name=d['name'])
