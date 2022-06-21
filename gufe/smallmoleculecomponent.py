# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import logging
# openff complains about oechem being missing, shhh
logger = logging.getLogger('openff.toolkit')
logger.setLevel(logging.ERROR)
from openff.toolkit.topology import Molecule as OFFMolecule
from openff.toolkit.utils.serialization import Serializable
import warnings

from rdkit import Chem

from gufe import __version__
from gufe import Component
from gufe.molhashing import hashmol
from gufe.custom_typing import RDKitMol, OEMol
from gufe.storage.utils import SerializationInfo


def _ensure_ofe_name(mol: RDKitMol, name: str) -> str:
    """
    Determine the correct name from the rdkit.Chem.Mol and the user-provided
    name; ensure that is set in the rdkit representation.
    """
    try:
        rdkit_name = mol.GetProp("_Name")
    except KeyError:
        rdkit_name = ""

    try:
        rdkit_name = mol.GetProp("ofe-name")
    except KeyError:
        pass

    if name and rdkit_name and rdkit_name != name:
        warnings.warn(f"SmallMoleculeComponent being renamed from {rdkit_name}"
                      f"to {name}.")
    elif name == "":
        name = rdkit_name

    mol.SetProp("ofe-name", name)
    return name


def _ensure_ofe_version(mol: RDKitMol):
    """Ensure the rdkit representation has the current version associated"""
    mol.SetProp("ofe-version", __version__)


class SmallMoleculeComponent(Component, Serializable):
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

    _storage_path = "setup/components/{md5}.sdf"

    def __init__(self, rdkit: RDKitMol, name: str = ""):
        name = _ensure_ofe_name(rdkit, name)
        _ensure_ofe_version(rdkit)
        self._rdkit = rdkit
        self._hash = hashmol(self._rdkit, name=name)

    def to_rdkit(self) -> RDKitMol:
        """Return an RDKit copied representation of this molecule"""
        return Chem.Mol(self._rdkit)

    @classmethod
    def from_rdkit(cls, rdkit: RDKitMol, name: str = ""):
        """Create a SmallMoleculeComponent, copying from an RDKit Mol"""
        return cls(rdkit=Chem.Mol(rdkit), name=name)

    def to_openeye(self) -> OEMol:
        """OEChem representation of this molecule"""
        return self.to_openff().to_openeye()

    @classmethod
    def from_openeye(cls, oemol: OEMol, name: str = ""):
        raise NotImplementedError

    def to_openff(self) -> OFFMolecule:
        """OpenFF Toolkit representation of this molecule"""
        m = OFFMolecule(self._rdkit, allow_undefined_stereo=True)
        m.name = self.name

        return m

    @classmethod
    def from_openff(cls, openff: OFFMolecule, name: str = ""):
        """Construct from an OpenFF toolkit Molecule"""
        return cls(openff.to_rdkit(), name=name)

    @property
    def smiles(self) -> str:
        return self._hash.smiles

    @property
    def name(self) -> str:
        return self._hash.name

    def __hash__(self):
        return hash(self._hash)

    def __eq__(self, other):
        return hash(self) == hash(other)

    def to_dict(self) -> dict:
        """Serialize to dict representation"""
        d = self.to_openff().to_dict()
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Deserialize from dict representation"""
        return cls.from_openff(OFFMolecule.from_dict(d),
                               name=d.get('name', ''))

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

    def to_storage_ready(self):
        bytes_data = self.to_sdf().encode("utf-8")
        metadata = SerializationInfo.create_metadata(self, bytes_data,
                                                     self._storage_path)
        return {self: SerializationInfo(bytes_data, metadata)}

    @classmethod
    def from_storage_bytes(cls, serialized_bytes, load_func):
        return cls.from_sdf_string(serialized_bytes.decode("utf-8"))

    @property
    def total_charge(self):
        return Chem.GetFormalCharge(self._rdkit)
