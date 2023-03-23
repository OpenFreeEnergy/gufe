import abc
import json
import warnings
from rdkit import Chem
from typing import Optional, Any

from .component import Component

# typing
from ..custom_typing import RDKitMol
from ..molhashing import serialize_numpy, deserialize_numpy


_INT_TO_ATOMCHIRAL = {
    0: Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
    1: Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
    2: Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW,
    3: Chem.rdchem.ChiralType.CHI_OTHER,
    4: Chem.rdchem.ChiralType.CHI_TETRAHEDRAL,
    5: Chem.rdchem.ChiralType.CHI_ALLENE,
    6: Chem.rdchem.ChiralType.CHI_SQUAREPLANAR,
    7: Chem.rdchem.ChiralType.CHI_TRIGONALBIPYRAMIDAL,
    8: Chem.rdchem.ChiralType.CHI_OCTAHEDRAL}
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
    """Convert from internal dict serialisation format to rdkit mol"""
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


class ExplicitMoleculeComponent(Component, abc.ABC):
    """Base class for explicit molecules.

    This provides basic serialization and conversion to different
    representations. Specific file formats, such as SDF for small molecules
    or PDB for proteins, should be implemented in subclasses.
    """
    _rdkit: Chem.Mol
    _name: str

    def __init__(self, rdkit: RDKitMol, name: str = ""):
        name = _ensure_ofe_name(rdkit, name)
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
