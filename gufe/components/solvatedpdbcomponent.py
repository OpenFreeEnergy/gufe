# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import warnings
from os import PathLike
from typing import TextIO

import numpy as np
from openff.toolkit import Quantity
from openff.units import unit as offunit
from openff.units.openmm import from_openmm
from openmm import unit as omm_unit
from rdkit import Chem
from rdkit.Chem.rdchem import Mol

from ..vendor.openff.interchange._annotations import _is_box_shape
from ..vendor.openff.interchange._packmol import _box_vectors_are_in_reduced_form
from ..vendor.pdb_file.pdbfile import PDBFile
from ..vendor.pdb_file.pdbxfile import PDBxFile
from .proteincomponent import ProteinComponent
from .solventcomponent import BaseSolventComponent


class SolvatedPDBComponent(ProteinComponent, BaseSolventComponent):
    """
    Protein component with explicit solvent and box vectors.

    This class represents a protein structure that is associated with
    explicit box vectors. Unlike ``ProteinComponent``, instances
    of this class always have box vectors, which are treated as
    part of the component's identity (affecting equality and hashing).

    Notes
    -----
    * ``box_vectors`` must be an OpenFF quantity with units.
    * Box vectors are serialized and included in equality and hash checks.
    * Construction will fail if box vectors cannot be determined as well as
      if the RDKit molecule has only one disconnected fragment.
    """

    def __init__(self, rdkit: Mol, box_vectors: Quantity, name: str = ""):
        """
        Parameters
        ----------
        rdkit : rdkit.Chem.Mol
            RDKit representation of the protein.
        box_vectors : openff.units.Quantity
            Periodic box vectors with units of length, compatible with
            nanometers. Must be a (3, 3) array in reduced form.
            Reduced form is a canonical representation of the unit cell and removes
            ambiguity in periodic boundary conditions (see
            https://docs.openmm.org/latest/userguide/theory/05_other_features.html).
        name : str, optional
            Name of the component.

        Raises
        ------
        TypeError
            If ``box_vectors`` is not an OpenFF Quantity.
        ValueError
            If ``box_vectors`` are not valid box vectors.
            If the RDKit molecule contains only one disconnected fragment.
        """
        self._validate_box_vectors(box_vectors)
        super().__init__(rdkit=rdkit, name=name)
        self._validate_multiple_molecules(rdkit)
        self.box_vectors = box_vectors
        # Density sanity check (warning only)
        self._warn_if_density_too_low()

    @staticmethod
    def _validate_box_vectors(box):
        """
        Validate box vectors.

        Parameters
        ----------
        box : openff.units.Quantity
            Box vectors to validate.

        Raises
        ------
        TypeError
            If ``box`` is not an OpenFF Quantity.
        ValueError
            If ``box`` does not represent valid reduced-form box vectors.
        """
        if box is None:
            raise ValueError("box_vectors must be provided")

        # OpenFF Quantity check
        _is_box_shape(box)

        # Reduced-form check
        if not _box_vectors_are_in_reduced_form(box):
            raise ValueError(f"box_vectors: {box} are not in OpenMM reduced form")

    @staticmethod
    def _is_water_fragment(mol: Mol) -> bool:
        """
        Return True if this fragment looks like a water molecule (TIP3P/TIP4P/etc).

        Definition:
        - exactly 1 oxygen
        - exactly 2 hydrogens
        Atoms with atomic number 0 (virtual sites) are ignored.
        """
        if mol.GetNumAtoms() != 3:
            return False
        n_H = 0
        n_O = 0

        for atom in mol.GetAtoms():
            match atom.GetAtomicNum():
                case 1:
                    n_H += 1
                case 8:
                    n_O += 1
                case _:
                    return False

        return (n_H == 2) and (n_O == 1)


    @classmethod
    def _count_waters(cls, rdkit_mol: Mol) -> int:
        """
        Count water molecules by disconnected fragments.
        """
        frags = Chem.rdmolops.GetMolFrags(rdkit_mol, asMols=True)

        return sum(cls._is_water_fragment(frag) for frag in frags)

    @classmethod
    def _validate_multiple_molecules(cls, rdkit_mol, *, min_waters: int = 50):
        """
        Ensure multiple fragments are present and warn if fewer than `min_waters`
        waters are detected.
        """
        frags = Chem.rdmolops.GetMolFrags(rdkit_mol, asMols=False)
        if len(frags) <= 1:
            raise ValueError(
                "SolvatedPDBComponent requires multiple molecules (e.g., protein + solvent). Found a single molecule."
            )

        n_waters = cls._count_waters(rdkit_mol)
        print(n_waters)

        if n_waters < min_waters:
            warnings.warn(
                f"Only {n_waters} water molecules detected (expected ≥ {min_waters}). "
                "This may indicate missing solvent or a non-aqueous system.",
                UserWarning,
            )

    def compute_density(self):
        """
        Estimate the system density in g/L from the RDKit molecule and box vectors.

        Returns
        -------
        density : openff.units.Quantity
            Estimated density in grams per liter.
        """
        # total mass
        total_mass = sum(atom.GetMass() for atom in self._rdkit.GetAtoms()) * offunit.dalton

        # box volume
        box_nm = self.box_vectors.to("nanometer").magnitude
        volume_nm3 = abs(np.linalg.det(box_nm)) * offunit.nanometer**3
        volume_L = volume_nm3.to("liter")

        # density
        density = total_mass.to("gram") / volume_L
        return density

    def _warn_if_density_too_low(self):
        """
        Give a warning if the estimated system density is below 500 g/l.
        This is a heuristic check intended to catch issues such as missing
        solvent or box vector that are too big.
        """
        density = self.compute_density()
        min_density = 500 * offunit.gram / offunit.liter

        if density < min_density:
            warnings.warn(
                "Estimated system density is very low.\n"
                f"  Density: {density:.3f}\n"
                "This usually indicates missing solvent or incorrect box "
                "vectors.\n",
                UserWarning,
            )

    @staticmethod
    def _estimate_box(omm_structure):
        """
        Estimate an orthorhombic box from atomic coordinates.
        The bounding box is computed from the minimum and maximum atomic
        coordinates and returned as orthorhombic box vectors.

        Parameters
        ----------
        omm_structure : PDBFile or PDBxFile
            OpenMM structure providing atomic positions.

        Returns
        -------
        openff.units.Quantity
            Orthorhombic box vectors with units of nanometers.
        """
        coords_nm = np.asarray(omm_structure.positions.value_in_unit(omm_unit.nanometer))

        mins = coords_nm.min(axis=0)
        maxs = coords_nm.max(axis=0)
        lengths = maxs - mins

        box = np.array(
            [
                [lengths[0], 0.0, 0.0],
                [0.0, lengths[1], 0.0],
                [0.0, 0.0, lengths[2]],
            ]
        )

        return box * offunit.nanometer

    @classmethod
    def _resolve_box_vectors(
        cls,
        structure,
        *,
        box_vectors=None,
        infer_box_vectors: bool = False,
    ):
        """
        Resolve periodic box vectors from user input, file, or inference.

        Parameters
        ----------
        structure : PDBFile or PDBxFile
            Loaded OpenMM structure.
        box_vectors : openff.units.Quantity, optional
            Explicit box vectors to associate with the component.
        infer_box_vectors : bool, optional
            If True, estimate box vectors when not present in file.

        Returns
        -------
        openff.units.Quantity
            Box vectors with units of nanometers.

        Raises
        ------
        ValueError
            If box vectors cannot be determined.
        """
        # 1. User-supplied box vectors win
        if box_vectors is not None:
            return box_vectors

        # 2. Try reading box vectors from the file
        box = structure.topology.getPeriodicBoxVectors()
        if box is not None:
            box = from_openmm(box)
            # Cryo-EM special case: unit cell present but meaningless (~1 Å³)
            # Treat this as "no box vectors"
            lengths = np.diag(box.to("nanometer").magnitude)
            if np.allclose(lengths, 0.1, atol=1e-3):
                warnings.warn(
                    "Periodic box vectors of ~1 Å detected in the input structure. "
                    "This can e.g. be observed in cryo-EM derived PDB files and does "
                    "not represent a meaningful simulation box. The box vectors will "
                    "be ignored.",
                    UserWarning,
                )
            else:
                return box

        # 3. Infer box vectors if requested
        if infer_box_vectors:
            box = cls._estimate_box(structure)
            warnings.warn(
                "Box vectors were inferred from the atomic coordinates.\n"
                "Note: This heuristic assumes that the coordinates reflect the true "
                "periodic unit cell. It may produce incorrect box dimensions for "
                "structures that were post-processed (e.g., unwrapped after MD).\n"
                f"Inferred box vectors:\n{box}",
                UserWarning,
            )
            return box

        raise ValueError(
            "Could not determine box_vectors. Please provide them explicitly "
            "via the ``box_vectors`` argument or enable ``infer_box_vectors``"
        )

    @classmethod
    def from_pdb_file(
        cls,
        pdb_file: PathLike | TextIO,
        name: str = "",
        *,
        box_vectors=None,
        infer_box_vectors: bool = False,
    ):
        """
        Create a SolvatedPDBComponent from a PDB file.
        """
        pdb = PDBFile(pdb_file)

        box = cls._resolve_box_vectors(
            pdb,
            box_vectors=box_vectors,
            infer_box_vectors=infer_box_vectors,
        )

        return cls._from_openmmPDBFile(
            pdb,
            name=name,
            box_vectors=box,
        )

    @classmethod
    def from_pdbx_file(
        cls,
        pdbx_file: str,
        name: str = "",
        *,
        box_vectors=None,
        infer_box_vectors: bool = False,
    ):
        """
        Create a SolvatedPDBComponent from a PDBx/mmCIF file.
        """
        pdbx = PDBxFile(pdbx_file)

        box = cls._resolve_box_vectors(
            pdbx,
            box_vectors=box_vectors,
            infer_box_vectors=infer_box_vectors,
        )

        return cls._from_openmmPDBFile(
            pdbx,
            name=name,
            box_vectors=box,
        )

    @classmethod
    def _from_openmmPDBFile(
        cls,
        openmm_PDBFile,
        *,
        name="",
        box_vectors,
    ):
        """
        Construct a SolvatedPDBComponent from an OpenMM PDBFile or PDBxFile.

        This method converts an OpenMM structure into the internal RDKit-based
        representation. Box vectors are required and must be present.

        Parameters
        ----------
        openmm_PDBFile : openmm.app.PDBFile or openmm.app.PDBxFile
            OpenMM object containing topology, positions, and (optionally)
            box vectors.
        name : str, optional
            Name of the protein component.
        box_vectors: openff.units.Quantity
            Box vectors with units of nanometers.

        Returns
        -------
        SolvatedPDBComponent
            The constructed solvated protein component.

        Raises
        ------
        ValueError
            If ``box_vectors`` are not provided.
        """
        if box_vectors is None:
            raise ValueError("Box vectors are required but were not provided.")

        prot = ProteinComponent._from_openmmPDBFile(openmm_PDBFile, name=name)

        return cls(
            rdkit=prot._rdkit,
            name=prot.name,
            box_vectors=box_vectors,
        )

    def _to_dict(self):
        """
        Serialize the component to a dictionary.

        Periodic box vectors are always serialized explicitly as a numeric
        array plus unit string.
        """
        d = super()._to_dict()
        box = self.box_vectors.to("nanometer")

        d["box_vectors"] = box

        return d

    @classmethod
    def _from_dict(cls, d, name=""):
        """
        Deserialize from a dictionary.
        """
        box_vectors = d.get("box_vectors")
        if box_vectors is None:
            raise ValueError("box_vectors must be present in the serialized dict")

        prot = ProteinComponent._from_dict(d.copy(), name=name)

        return cls(
            rdkit=prot._rdkit,
            name=prot.name,
            box_vectors=box_vectors,
        )


class ProteinMembraneComponent(SolvatedPDBComponent):
    """
    Solvated protein component with an explicit membrane and periodic box.

    This class inherits all behavior and
    constraints from ``SolvatedPDBComponent`` and does not introduce any
    additional data or methods.
    The primary purpose of this subclass is semantic: it acts as a marker
    indicating the presence of a membrane. Code elsewhere may use this
    distinction (e.g., via ``isinstance`` checks or type annotations) to
    enable membrane-specific behavior such as selecting a membrane-aware
    barostat or simulation protocol.

    Notes
    -----
    * This class inherits all behavior and requirements from
      ``SolvatedPDBComponent``.
    * This class does not add new fields or override behavior.
    * All requirements and guarantees of ``SolvatedPDBComponent`` apply.
    * The distinction between membrane and non-membrane systems is
      conveyed solely through the component type.
    """

    ...
