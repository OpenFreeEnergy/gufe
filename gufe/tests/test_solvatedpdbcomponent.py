# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe
import gzip
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_almost_equal
from openff.units import unit as offunit
from rdkit import Chem

from gufe import ProteinComponent, ProteinMembraneComponent, SolvatedPDBComponent

from ..vendor.openff.interchange._annotations import _is_box_shape
from ..vendor.openff.interchange._packmol import _box_vectors_are_in_reduced_form
from .conftest import OPENMM_VERSION
from .test_explicitmoleculecomponent import ExplicitMoleculeComponentMixin
from .test_tokenization import GufeTokenizableTestsMixin

if OPENMM_VERSION:
    from openmm import unit
    from openmm.app import pdbfile, pdbxfile


class TestSolvatedPDBComponent(GufeTokenizableTestsMixin, ExplicitMoleculeComponentMixin):
    cls = SolvatedPDBComponent
    repr = "SolvatedPDBComponent(name=Steve)"

    @pytest.fixture(scope="session")
    def instance(self, PDB_a2a_path):
        with gzip.open(PDB_a2a_path, "rb") as gzf:
            yield self.cls.from_pdb_file(gzf, name="Steve")

    def test_from_pdb_file_sets_box_vectors(self, instance):
        box = instance.box_vectors
        _is_box_shape(box)
        assert _box_vectors_are_in_reduced_form(box)
        assert box[0, 0].m_as(offunit.nanometer) == pytest.approx(6.9587)

    def test_requires_box_vectors(self, PDB_a2a_path):
        with gzip.open(PDB_a2a_path, "rb") as gzf:
            prot = ProteinComponent.from_pdb_file(gzf)

        with pytest.raises(ValueError, match="box_vectors must be provided"):
            self.cls(
                rdkit=prot._rdkit,
                name=prot.name,
                box_vectors=None,
            )

    @pytest.mark.parametrize(
        "factory,path_fixture",
        [
            (SolvatedPDBComponent.from_pdb_file, "PDB_181L_path"),
            (SolvatedPDBComponent.from_pdbx_file, "PDBx_181L_path"),
        ],
    )
    def test_file_without_box_vectors_raises(self, factory, path_fixture, request):
        path = request.getfixturevalue(path_fixture)
        with pytest.raises(ValueError, match="Could not determine box_vectors"):
            factory(path)

    def test_box_vectors_preserved_in_dict_roundtrip(self, instance):
        d = instance.to_dict()
        m2 = self.cls.from_dict(d)

        v1 = instance.box_vectors.m
        v2 = m2.box_vectors.m

        assert_almost_equal(actual=v1, desired=v2, decimal=6)

    def test_from_dict_requires_box_vectors(self, instance):
        d = instance.to_dict()

        # Remove box vectors to simulate invalid / legacy data
        d.pop("box_vectors")

        with pytest.raises(
            ValueError,
            match="box_vectors must be present in the serialized dict",
        ):
            instance.__class__.from_dict(d)

    @pytest.mark.parametrize(
        "factory,loader,path_fixture",
        [
            (
                SolvatedPDBComponent.from_pdb_file,
                pdbfile.PDBFile,
                "PDB_a2a_path",
            ),
            (
                SolvatedPDBComponent.from_pdbx_file,
                pdbxfile.PDBxFile,
                "PDBx_a2a_path",
            ),
        ],
    )
    def test_uses_file_box_vectors(self, factory, loader, path_fixture, request):
        path = request.getfixturevalue(path_fixture)

        with gzip.open(path, "rt") as f:
            comp = factory(f)
        with gzip.open(path, "rt") as f2:
            ref = loader(f2).topology.getPeriodicBoxVectors()

        actual = comp.box_vectors.m_as(offunit.nanometer)
        expected = ref.value_in_unit(unit.nanometer)

        assert_almost_equal(actual=actual, desired=expected, decimal=6)

    @pytest.mark.parametrize(
        "factory,path_fixture",
        [
            (SolvatedPDBComponent.from_pdb_file, "PDB_181L_path"),
            (SolvatedPDBComponent.from_pdbx_file, "PDBx_181L_path"),
        ],
    )
    def test_infer_box_vectors_produces_valid_box(self, factory, path_fixture, request):
        path = request.getfixturevalue(path_fixture)

        with pytest.warns(UserWarning, match="Box vectors were inferred"):
            comp = factory(path, infer_box_vectors=True)
            box = comp.box_vectors

            _is_box_shape(box)
            assert _box_vectors_are_in_reduced_form(box)

    @pytest.mark.parametrize(
        "factory,path_fixture",
        [
            (SolvatedPDBComponent.from_pdb_file, "PDB_a2a_path"),
            (SolvatedPDBComponent.from_pdbx_file, "PDBx_a2a_path"),
        ],
    )
    def test_realistic_system_density(self, factory, path_fixture, request):
        """
        Check that a properly solvated system produces a realistic density.
        """
        path = request.getfixturevalue(path_fixture)
        with gzip.open(path, "rt") as f:
            comp = factory(f)

        density = comp.density
        # Expect realistic protein + solvent density ~> 800-1200 g/L (0.8-1.2 g/mL)
        assert density > 800 * offunit.gram / offunit.liter
        assert density < 1200 * offunit.gram / offunit.liter

    @pytest.mark.parametrize(
        "factory,path_fixture",
        [
            (SolvatedPDBComponent.from_pdb_file, "PDB_a2a_path"),
            (SolvatedPDBComponent.from_pdbx_file, "PDBx_a2a_path"),
        ],
    )
    def test_validate_passes(self, factory, path_fixture, request):
        path = request.getfixturevalue(path_fixture)
        with gzip.open(path, "rt") as f:
            comp = factory(f)
        # Should not raise for proper box vector
        comp.validate()

    def test_validate_low_density_raises(self, instance):
        # Inflate the box to reduce density
        assert instance.density > 500 * offunit.gram / offunit.liter
        instance.box_vectors = instance.box_vectors * 100
        assert instance.density < 500 * offunit.gram / offunit.liter
        with pytest.raises(ValueError, match="density is very low"):
            instance.validate(min_density=500 * offunit.gram / offunit.liter)

    def test_box_vectors_affect_equality(self, instance):
        v = np.eye(3) * 2.0 * offunit.nanometer

        comp2 = instance.copy_with_replacements(box_vectors=v)

        assert instance != comp2

    def test_from_openmmPDBFile_raises_without_box_vectors(self, PDB_181L_path):
        pdb = pdbfile.PDBFile(PDB_181L_path)

        with pytest.raises(ValueError, match="Box vectors are required"):
            self.cls._from_openmmPDBFile(
                pdb,
                name="test",
                box_vectors=None,
            )

    def test_from_pdbx_file_user_box_vectors(self, PDBx_a2a_path):
        b = np.eye(3) * 2.0 * offunit.nanometer
        with gzip.open(PDBx_a2a_path, "rt") as f:
            comp = self.cls.from_pdbx_file(f, box_vectors=b)
        box = comp.box_vectors
        assert box is not None
        assert box.shape == (3, 3)
        assert box.units.is_compatible_with(offunit.nanometer)

    @pytest.mark.parametrize(
        "factory,path_fixture",
        [
            (SolvatedPDBComponent.from_pdb_file, "PDB_a2a_path"),
            (SolvatedPDBComponent.from_pdbx_file, "PDBx_a2a_path"),
        ],
    )
    def test_explicit_box_vectors_override_file_box(self, factory, path_fixture, request):
        path = request.getfixturevalue(path_fixture)

        with gzip.open(path, "rt") as f:
            ref = factory(f)
        with gzip.open(path, "rt") as f2:
            override = np.eye(3) * 2.0 * offunit.nanometer
            comp = factory(f2, box_vectors=override)
        ref_box = ref.box_vectors

        assert not np.allclose(
            ref_box.m_as(offunit.nanometer),
            override.m_as(offunit.nanometer),
        )

        assert_almost_equal(
            comp.box_vectors.m_as(offunit.nanometer),
            override.m_as(offunit.nanometer),
            decimal=6,
        )

    def test_box_vectors_not_reduced_form(self, instance):
        bad = (
            np.array(
                [
                    [2.0, 0.0, 0.0],
                    [3.0, 2.0, 0.0],  # invalid reduced form
                    [0.0, 0.0, 2.0],
                ]
            )
            * offunit.nanometer
        )

        with pytest.raises(ValueError, match="reduced form"):
            instance.copy_with_replacements(box_vectors=bad)

    def test_cryo_em_dummy_box(self, PDB_181L_path, tmp_path):
        pdb_text = Path(PDB_181L_path).read_text()

        # Insert a CRYST1 record with a 1 Å box
        cryo_em_cryst1 = "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n"
        pdb_path = tmp_path / "cryo_em.pdb"
        pdb_path.write_text(cryo_em_cryst1 + pdb_text)

        pdb = pdbfile.PDBFile(str(pdb_path))

        with pytest.warns(UserWarning, match="cryo-EM"):
            with pytest.raises(ValueError, match="box_vectors"):
                SolvatedPDBComponent._resolve_box_vectors(pdb)


class TestProteinMembraneComponent(GufeTokenizableTestsMixin, ExplicitMoleculeComponentMixin):
    cls = ProteinMembraneComponent
    repr = "ProteinMembraneComponent(name=Steve)"

    @pytest.fixture(scope="session")
    def instance(self, PDB_a2a_path):
        with gzip.open(PDB_a2a_path, "rb") as gzf:
            yield self.cls.from_pdb_file(gzf, name="Steve")

    @pytest.mark.parametrize(
        "factory,path_fixture",
        [
            (ProteinMembraneComponent.from_pdb_file, "PDB_a2a_path"),
            (ProteinMembraneComponent.from_pdbx_file, "PDBx_a2a_path"),
        ],
    )
    def test_validate_passes(self, factory, path_fixture, request):
        path = request.getfixturevalue(path_fixture)
        with gzip.open(path, "rt") as f:
            comp = factory(f)
        # Should not raise for properly solvated system
        comp.validate(min_waters=50)

    def test_is_water_fragment_unknown_atom(self):
        # Build a fragment with 3 atoms: 1 O, 1 H, 1 C (C triggers `case _`)
        mol = Chem.RWMol()
        o = mol.AddAtom(Chem.Atom(8))  # oxygen
        h = mol.AddAtom(Chem.Atom(1))  # hydrogen
        c = mol.AddAtom(Chem.Atom(6))  # carbon

        mol.AddBond(o, h, Chem.BondType.SINGLE)
        mol.AddBond(o, c, Chem.BondType.SINGLE)

        mol = mol.GetMol()

        assert ProteinMembraneComponent._is_water_fragment(mol) is False

    def test_validate_few_waters_raises(self, PDB_181L_path):
        comp = ProteinMembraneComponent.from_pdb_file(PDB_181L_path, infer_box_vectors=True)
        # Only 8 Xray waters, not properly solvated
        assert comp.n_waters == 8
        with pytest.raises(ValueError, match="water molecules detected"):
            comp.validate(min_waters=50)

    def test_validate_few_waters_and_low_density_raises(self, PDB_181L_path):
        comp = ProteinMembraneComponent.from_pdb_file(PDB_181L_path, infer_box_vectors=True)

        with pytest.raises(ValueError) as excinfo:
            comp.validate()

        # Check that both error messages appear
        err_msg = str(excinfo.value)
        assert "water molecules detected" in err_msg
        assert "density is very low" in err_msg
