import json
import warnings
from rdkit import Chem

from .component import Component
from .. import __version__
from ..molhashing import hashmol

# typing
from ..custom_typing import RDKitMol


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


class ExplicitMoleculeComponent(Component):
    """Base class for explicit molecules.

    This provides basic serialization and conversion to different
    representations. Specific file formats, such as SDF for small molecules
    or PDB for proteins, should be implemented in subclasses.

    We default to use SMILES as the basis of the hash, but this can be
    overridden using the `_hashmol` method.
    """
    def __init__(self, rdkit: RDKitMol, name: str = ""):
        name = _ensure_ofe_name(rdkit, name)
        _ensure_ofe_version(rdkit)
        conformers = list(rdkit.GetConformers())
        if not conformers:
            raise ValueError("Molecule was provided with no conformers.")
        
        n_confs = len(conformers)
        if n_confs > 1:
            warnings.warn(f"Molecule provided with {n_confs} conformers. "
                          "Only the first will be used.")

        self._rdkit = rdkit
        self._hash = self._hashmol(name=name)

    def __hash__(self):
        return hash(self._hash)

    def __eq__(self, other):
        return hash(self) == hash(other)

    def _hashmol(self, name):
        return hashmol(self._rdkit, name=name)

    def _defaults(self):
        return super()._defaults()

    # Props
    @property
    def name(self) -> str:
        return self._hash.name

    @property
    def smiles(self) -> str:
        return self._hash.smiles

    @property
    def total_charge(self):
        return Chem.GetFormalCharge(self._rdkit)

    # TO
    def to_rdkit(self) -> RDKitMol:
        """Return an RDKit copied representation of this molecule"""
        return Chem.Mol(self._rdkit)

    def to_json(self):
        return json.dumps(self.to_dict())

    # From
    @classmethod
    def from_json(cls, json_str):
        dct = json.loads(json_str)
        return cls.from_dict(dct)

    @classmethod
    def from_rdkit(cls, rdkit: RDKitMol, name: str = ""):
        """Create a SmallMoleculeComponent, copying from an RDKit Mol"""
        return cls(rdkit=Chem.Mol(rdkit), name=name)

    # Not Implemented - interface Functions:
    @classmethod
    def _from_dict(cls, d: dict):
        """Deserialize from dict representation"""
        raise NotImplemented()

    def _to_dict(self) -> dict:
        raise NotImplemented()
