import json
import numpy as np
import warnings
from rdkit import Chem
from typing import Any, Optional

from .component import Component

from ..custom_typing import RDKitMol
from ..molhashing import deserialize_numpy, serialize_numpy


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
        warnings.warn(
            f"Component being renamed from {rdkit_name}"
            f"to {name}."
        )
    elif name == "":
        name = rdkit_name

    mol.SetProp("ofe-name", name)
    return name


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


def _mol_from_dict(d: dict) -> Chem.Mol:
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

    return m


def _mol_to_dict(m: Chem.Mol) -> dict[str, Any]:
    # in a perfect world we'd use ToBinary()
    # but this format slowly evolves, so the future hash of a SMC could change if rdkit were updated
    # this is based on that method, with some irrelevant fields cut out

    output: dict[str, Any] = {}

    atoms = []
    for atom in m.GetAtoms():
        atoms.append((
            atom.GetAtomicNum(), atom.GetIsotope(), atom.GetFormalCharge(), atom.GetIsAromatic(),
            _ATOMCHIRAL_TO_INT[atom.GetChiralTag()], atom.GetAtomMapNum(),
            atom.GetPropsAsDict(includePrivate=False),
        ))
    output['atoms'] = atoms

    bonds = []
    for bond in m.GetBonds():
        bonds.append((
            bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), _BONDTYPE_TO_INT[bond.GetBondType()],
            _BONDSTEREO_TO_INT[bond.GetStereo()],
            bond.GetPropsAsDict(includePrivate=False)
        ))
    output['bonds'] = bonds

    conf = m.GetConformer()
    output['conformer'] = (serialize_numpy(conf.GetPositions()), conf.GetPropsAsDict(includePrivate=False))

    output['molprops'] = m.GetPropsAsDict(includePrivate=False)

    return output


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
    if 'atom.dprop.PartialCharge' not in mol.GetPropNames():
        return

    p_chgs = np.array(
        mol.GetProp('atom.dprop.PartialCharge').split(), dtype=float
    )

    if len(p_chgs) != mol.GetNumAtoms():
        errmsg = (f"Incorrect number of partial charges: {len(p_chgs)} "
                  f" were provided for {mol.GetNumAtoms()} atoms")
        raise ValueError(errmsg)

    if (sum(p_chgs) - Chem.GetFormalCharge(mol)) > 0.01:
        errmsg = (f"Sum of partial charges {sum(p_chgs)} differs from "
                  f"RDKit formal charge {Chem.GetFormalCharge(mol)}")
        raise ValueError(errmsg)

    if np.all(np.isclose(p_chgs, 0.0)):
        wmsg = (f"Partial charges provided all equal to "
                "zero. These may be ignored by some Protocols.")
        warnings.warn(wmsg)
    else:
        wmsg = ("Partial charges have been provided, these will "
                "preferentially be used instead of generating new "
                "partial charges")
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
            warnings.warn(
                f"Molecule provided with {n_confs} conformers. "
                f"Only the first will be used."
            )

        if not any(atom.GetAtomicNum() == 1 for atom in rdkit.GetAtoms()):
            warnings.warn("Molecule doesn't have any hydrogen atoms present. "
                          "If this is unexpected, consider loading the molecule with `removeHs=False`")

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

    def to_rdkit(self) -> RDKitMol:
        """Return an RDKit copied representation of this molecule"""
        return Chem.Mol(self._rdkit)

    def to_json(self):
        return json.dumps(self.to_dict())

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
