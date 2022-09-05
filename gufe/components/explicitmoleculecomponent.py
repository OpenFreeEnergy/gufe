from .component import Component
import warnings
from rdkit import Chem
from .. import __version__
from ..molhashing import hashmol, deserialize_numpy, serialize_numpy
from openff.units import unit

# typing
from openff.toolkit.topology import Molecule as OFFMolecule
from ..custom_typing import RDKitMol, OEMol


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
        elif (n_confs := len(conformers)) > 1:
            warnings.warn(f"Molecule provided with {n_confs} conformers. "
                          "Only the first will be used.")

        self._rdkit = rdkit
        self._hash = self._hashmol(name=name)

    def _hashmol(self, name):
        return hashmol(self._rdkit, name=name)

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
        # required attributes: (based on openff to_dict)
        # for each atom:
        #   element, name, formal charge, aromaticity, stereochemistry
        # for each bond:
        #   idx0, idx1, order, aromaticity, stereochemistry
        # TODO: Do we care about fractional bond orders?
        #       is aromaticity reperceived on creation?
        # NOTE: Here we're implicitly using units of angstrom and elementary
        # charge. We might want to explcitly include them in the stored dict.
        m = self.to_openff()
        atoms = [
            (atom.atomic_number,
             atom.name,
             atom.formal_charge.m_as(unit.elementary_charge),
             atom.is_aromatic,
             atom.stereochemistry or '')
            for atom in m.atoms
        ]
        bonds = [
            (bond.atom1_index, bond.atom2_index, bond.bond_order,
             bond.is_aromatic, bond.stereochemistry or '')
            for bond in m.bonds
        ]

        if m.conformers is None:  # -no-cov-
            # this should not be reachable; indicates that something went
            # very wrong
            raise RuntimeError(f"{self.__class__.__name__} must have at "
                               "least 1 conformer")

        conformers = [
            serialize_numpy(conf.m_as(unit.angstrom))
            for conf in m.conformers
        ]

        d = {
            'atoms': atoms,
            'bonds': bonds,
            'name': self.name,
            'conformers': conformers,
        }

        return d

    @classmethod
    def _from_dict(cls, d: dict):
        """Deserialize from dict representation"""
        # manually construct OpenFF molecule as in cookbook
        m = OFFMolecule()
        for (an, name, fc, arom, stereo) in d['atoms']:
            m.add_atom(
                atomic_number=an,
                formal_charge=fc * unit.elementary_charge,
                is_aromatic=arom,
                stereochemistry=stereo or None,
                name=name,
            )

        for (idx1, idx2, order, arom, stereo) in d['bonds']:
            m.add_bond(
                atom1=idx1,
                atom2=idx2,
                bond_order=order,
                is_aromatic=arom,
                stereochemistry=stereo or None,
            )

        for conf in d['conformers']:
            m.add_conformer(deserialize_numpy(conf) * unit.angstrom)

        return cls.from_openff(m, name=d['name'])

    def _defaults(self):
        return super()._defaults()

    def __hash__(self):
        return hash(self._hash)

    def __eq__(self, other):
        return hash(self) == hash(other)

    @property
    def name(self) -> str:
        return self._hash.name

    @property
    def smiles(self) -> str:
        return self._hash.smiles

    def to_json(self):
        return json.dumps(self.to_dict())

    @classmethod
    def from_json(cls, json_str):
        dct = json.loads(json_str)
        return cls.from_dict(dct)

    @property
    def total_charge(self):
        return Chem.GetFormalCharge(self._rdkit)
