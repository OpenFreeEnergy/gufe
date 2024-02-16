import json
import warnings
from typing import Optional

import numpy as np
from rdkit import Chem

# typing
from ..custom_typing import RDKitMol
from .component import Component


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
        warnings.warn(f"Component being renamed from {rdkit_name}" f"to {name}.")
    elif name == "":
        name = rdkit_name

    mol.SetProp("ofe-name", name)
    return name


def _check_partial_charges(mol: RDKitMol) -> None:
    """
    Checks for the presence of partial charges.

    Raises
    ------
    ValueError
      * If the partial charges are not of length atoms.
      * If the sum of partial charges is not equal to the
        formal charge.
    UserWarning
      * If partial charges are found.
      * If the partial charges are near 0 for all atoms.
    """
    if "atom.dprop.PartialCharge" not in mol.GetPropNames():
        return

    p_chgs = np.array(mol.GetProp("atom.dprop.PartialCharge").split(), dtype=float)

    if len(p_chgs) != mol.GetNumAtoms():
        errmsg = f"Incorrect number of partial charges: {len(p_chgs)} " f" were provided for {mol.GetNumAtoms()} atoms"
        raise ValueError(errmsg)

    if (sum(p_chgs) - Chem.GetFormalCharge(mol)) > 0.01:
        errmsg = (
            f"Sum of partial charges {sum(p_chgs)} differs from " f"RDKit formal charge {Chem.GetFormalCharge(mol)}"
        )
        raise ValueError(errmsg)

    if np.all(np.isclose(p_chgs, 0.0)):
        wmsg = f"Partial charges provided all equal to " "zero. These may be ignored by some Protocols."
        warnings.warn(wmsg)
    else:
        wmsg = (
            "Partial charges have been provided, these will "
            "preferentially be used instead of generating new "
            "partial charges"
        )
        warnings.warn(wmsg)


class ExplicitMoleculeComponent(Component):
    """Base class for explicit molecules.

    This provides basic serialization and conversion to different
    representations. Specific file formats, such as SDF for small molecules
    or PDB for proteins, should be implemented in subclasses.
    """

    _rdkit: Chem.Mol
    _name: str

    def __init__(self, rdkit: RDKitMol, name: str = ""):
        name = _ensure_ofe_name(rdkit, name)
        _check_partial_charges(rdkit)
        conformers = list(rdkit.GetConformers())
        if not conformers:
            raise ValueError("Molecule was provided with no conformers.")

        n_confs = len(conformers)
        if n_confs > 1:
            warnings.warn(f"Molecule provided with {n_confs} conformers. " f"Only the first will be used.")

        if not any(atom.GetAtomicNum() == 1 for atom in rdkit.GetAtoms()):
            warnings.warn(
                "Molecule doesn't have any hydrogen atoms present. "
                "If this is unexpected, consider loading the molecule with `removeHs=False`",
            )

        self._rdkit = rdkit
        self._smiles: Optional[str] = None
        self._name = name

    @classmethod
    def _defaults(cls):
        return super()._defaults()

    # Props
    @property
    def name(self) -> str:
        return self._name

    @property
    def smiles(self) -> str:
        if self._smiles is None:
            self._smiles = Chem.MolToSmiles(Chem.RemoveHs(self._rdkit))

        return self._smiles

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
        """Create a Component, copying from an RDKit Mol"""
        return cls(rdkit=Chem.Mol(rdkit), name=name)

    # Not Implemented - interface Functions:
    @classmethod
    def _from_dict(cls, d: dict):
        """Deserialize from dict representation"""
        raise NotImplementedError()

    def _to_dict(self) -> dict:
        raise NotImplementedError()
