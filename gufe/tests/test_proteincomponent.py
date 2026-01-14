# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe
import copy
import io
import os
from pathlib import Path
from unittest import mock

import numpy as np
import pytest
from numpy.testing import assert_almost_equal
from openff.units import unit as offunit
from packaging.version import Version
from rdkit import Chem

from gufe import ProteinComponent, ProteinMembraneComponent, SolvatedPDBComponent

from ..molhashing import serialize_numpy
from .conftest import ALL_PDB_LOADERS, OPENMM_VERSION
from ..vendor.openff.interchange._annotations import _is_box_shape
from ..vendor.openff.interchange._packmol import _box_vectors_are_in_reduced_form
from .test_explicitmoleculecomponent import ExplicitMoleculeComponentMixin
from .test_tokenization import GufeTokenizableTestsMixin

if OPENMM_VERSION:
    from openmm import unit
    from openmm.app import pdbfile, pdbxfile


@pytest.fixture
def PDB_181L_mutant(PDB_181L_path):
    # has a single resname flipped to change the resulting sequence
    rdm = Chem.MolFromPDBFile(PDB_181L_path, removeHs=False)

    # important to select the Ca atom - Ca of residue 2
    at = rdm.GetAtomWithIdx(7)
    at.GetMonomerInfo().SetResidueName("X")

    return ProteinComponent.from_rdkit(rdm)


@pytest.fixture
def custom_pdb_ion(PDB_181L_path):
    def _make_custom_pdb_ion(new_ion: str):
        with open(PDB_181L_path, "r") as f:
            orig_pdb = f.read()

        str_to_replace = "HETATM 2615 CL  "  #  CL S 173      43.141  16.447   1.769  1.00  0.00          CL"

        new_str = str_to_replace.replace("CL  ", new_ion)

        test_pdb = orig_pdb.replace(str_to_replace, new_str)

        with io.StringIO(test_pdb) as f:
            yield f

    yield _make_custom_pdb_ion


def assert_same_pdb_lines(in_file_path, out_file_path):
    in_lines = []
    if hasattr(in_file_path, "readlines"):
        in_file = in_file_path
    else:
        in_file = open(in_file_path)

    if isinstance(out_file_path, io.StringIO):
        out_file = out_file_path
        must_close = False
    else:
        out_file = open(out_file_path)
        must_close = True

    in_lines = in_file.readlines()
    out_lines = out_file.readlines()

    if must_close:
        out_file.close()

    in_lines = [l for l in in_lines if not l.startswith(("REMARK", "CRYST", "# Created with"))]
    out_lines = [l for l in out_lines if not l.startswith(("REMARK", "CRYST", "# Created with"))]

    assert in_lines == out_lines


def assert_topology_equal(ref_top, top):
    assert ref_top.getNumAtoms() == top.getNumAtoms()
    for ref_atom, atom in zip(ref_top.atoms(), top.atoms()):
        assert ref_atom.name == atom.name
        assert ref_atom.element == atom.element
        assert ref_atom.index == atom.index
        assert ref_atom.id == atom.id
        assert ref_atom.residue.id == atom.residue.id

    assert ref_top.getNumResidues() == top.getNumResidues()
    for ref_res, res in zip(ref_top.residues(), top.residues()):
        assert len(ref_res) == len(res)
        assert ref_res.name == res.name
        assert ref_res.id == res.id
        assert ref_res.insertionCode == res.insertionCode
        assert ref_res.chain.id == res.chain.id

    for ref_chain, chain in zip(ref_top.chains(), top.chains()):
        assert len(ref_chain) == len(chain)
        assert ref_chain.id == chain.id

    assert ref_top.getNumBonds() == top.getNumBonds()
    for ref_bond, bond in zip(ref_top.bonds(), top.bonds()):
        assert ref_bond[0].index == bond[0].index
        assert ref_bond[1].index == bond[1].index


class TestProteinComponent(GufeTokenizableTestsMixin, ExplicitMoleculeComponentMixin):
    cls = ProteinComponent
    repr = "ProteinComponent(name=Steve)"

    @pytest.fixture(scope="session")
    def instance(self, PDB_181L_path):
        return self.cls.from_pdb_file(PDB_181L_path, name="Steve")

    # From
    @pytest.mark.parametrize("in_pdb_path", ALL_PDB_LOADERS.keys())
    def test_from_pdb_file(self, in_pdb_path):
        in_pdb_io = ALL_PDB_LOADERS[in_pdb_path]()
        p = self.cls.from_pdb_file(in_pdb_io, name="Steve")

        assert isinstance(p, ProteinComponent)
        assert p.name == "Steve"

    def test_from_pdbx_file(self, PDBx_181L_path):
        p = self.cls.from_pdbx_file(str(PDBx_181L_path), name="Steve")

        assert isinstance(p, ProteinComponent)
        assert p.name == "Steve"
        assert p.to_rdkit().GetNumAtoms() == 2639

    def test_from_rdkit(self, PDB_181L_path):
        m = Chem.MolFromPDBFile(PDB_181L_path, removeHs=False)
        p = self.cls.from_rdkit(rdkit=m, name="Steve")

        assert isinstance(p, ProteinComponent)
        assert p.name == "Steve"
        assert p.to_rdkit().GetNumAtoms() == 2639

    # To
    def test_to_rdkit(self, PDB_181L_path):
        pm = self.cls.from_pdb_file(PDB_181L_path)
        rdkitmol = pm.to_rdkit()

        assert isinstance(rdkitmol, Chem.Mol)
        assert rdkitmol.GetNumAtoms() == 2639

    def _test_file_output(self, input_path, output_path, input_type, output_func):
        if input_type == "filename":
            inp = str(output_path)
        elif input_type == "Path":
            inp = output_path
        elif input_type == "StringIO":
            inp = io.StringIO()
        elif input_type == "TextIOWrapper":
            inp = open(output_path, mode="w")

        output_func(inp)

        if input_type == "StringIO":
            inp.seek(0)
            output_path = inp

        if input_type == "TextIOWrapper":
            inp.close()

        assert_same_pdb_lines(in_file_path=str(input_path), out_file_path=output_path)

    @pytest.mark.parametrize("input_type", ["filename", "Path", "StringIO", "TextIOWrapper"])
    def test_to_pdbx_file(self, PDBx_181L_openMMClean_path, tmp_path, input_type):
        p = self.cls.from_pdbx_file(str(PDBx_181L_openMMClean_path), name="Bob")
        out_file_name = "tmp_181L_pdbx.cif"
        out_file = tmp_path / out_file_name

        self._test_file_output(
            input_path=PDBx_181L_openMMClean_path,
            output_path=out_file,
            input_type=input_type,
            output_func=p.to_pdbx_file,
        )

    @pytest.mark.parametrize("input_type", ["filename", "Path", "StringIO", "TextIOWrapper"])
    @pytest.mark.skipif(OPENMM_VERSION < Version("8.2"), reason="OpenMM version too old")
    @pytest.mark.skipif(OPENMM_VERSION == Version("0.0.0"), reason="OpenMM not installed")
    def test_to_pdb_input_types(self, PDB_181L_path, tmp_path, input_type):
        p = self.cls.from_pdb_file(str(PDB_181L_path), name="Bob")

        self._test_file_output(
            input_path=PDB_181L_path,
            output_path=tmp_path / "tmp_181L.pdb",
            input_type=input_type,
            output_func=p.to_pdb_file,
        )

    @pytest.mark.parametrize("in_pdb_path", ALL_PDB_LOADERS.keys())
    @pytest.mark.skipif(OPENMM_VERSION == Version("0.0.0"), reason="OpenMM not installed")
    def test_to_pdb_round_trip(self, in_pdb_path, tmp_path):
        in_pdb_io = ALL_PDB_LOADERS[in_pdb_path]()

        p = self.cls.from_pdb_file(in_pdb_io, name="Wuff")
        out_file_name = "tmp_" + in_pdb_path + ".pdb"
        out_file = tmp_path / out_file_name

        p.to_pdb_file(str(out_file))

        ref_in_pdb_io = ALL_PDB_LOADERS[in_pdb_path]()

        # generate openMM reference file:
        openmm_pdb = pdbfile.PDBFile(ref_in_pdb_io)
        out_ref_file_name = "tmp_" + in_pdb_path + "_openmm_ref.pdb"
        out_ref_file = tmp_path / out_ref_file_name

        pdbfile.PDBFile.writeFile(openmm_pdb.topology, openmm_pdb.positions, file=open(str(out_ref_file), "w"))
        assert_same_pdb_lines(in_file_path=str(out_ref_file), out_file_path=out_file)

    @pytest.mark.skipif(OPENMM_VERSION < Version("8.2"), reason="OpenMM version too old")
    def test_io_pdb_comparison(self, PDB_181L_path, tmp_path):
        out_file_name = "tmp_" + os.path.basename(PDB_181L_path)
        out_file = tmp_path / out_file_name

        p = self.cls.from_pdb_file(PDB_181L_path, name="Bob")
        _ = p.to_pdb_file(str(out_file))

        assert_same_pdb_lines(PDB_181L_path, str(out_file))

    def test_dummy_from_dict(self, PDB_181L_path):
        p = self.cls.from_pdb_file(PDB_181L_path, name="Bob")
        gufe_dict = p.to_dict()
        p2 = self.cls.from_dict(gufe_dict)

        assert p == p2

    # parametrize
    @pytest.mark.parametrize("in_pdb_path", ALL_PDB_LOADERS.keys())
    @pytest.mark.skipif(not OPENMM_VERSION, reason="OpenMM not installed")
    def test_to_openmm_positions(self, in_pdb_path):
        in_pdb_io = ALL_PDB_LOADERS[in_pdb_path]()
        ref_in_pdb_io = ALL_PDB_LOADERS[in_pdb_path]()

        openmm_pdb = pdbfile.PDBFile(ref_in_pdb_io)
        openmm_pos = openmm_pdb.positions

        p = self.cls.from_pdb_file(in_pdb_io, name="Bob")
        gufe_openmm_pos = p.to_openmm_positions()

        v1 = gufe_openmm_pos.value_in_unit(unit.nanometer)
        v2 = openmm_pos.value_in_unit(unit.nanometer)

        assert_almost_equal(actual=v1, desired=v2, decimal=6)

    # parametrize
    @pytest.mark.parametrize("in_pdb_path", ALL_PDB_LOADERS.keys())
    @pytest.mark.skipif(not OPENMM_VERSION, reason="OpenMM not installed")
    def test_to_openmm_topology(self, in_pdb_path):
        in_pdb_io = ALL_PDB_LOADERS[in_pdb_path]()
        ref_in_pdb_io = ALL_PDB_LOADERS[in_pdb_path]()

        openmm_pdb = pdbfile.PDBFile(ref_in_pdb_io)
        openmm_top = openmm_pdb.topology

        p = self.cls.from_pdb_file(in_pdb_io, name="Bob")
        gufe_openmm_top = p.to_openmm_topology()
        assert_topology_equal(openmm_top, gufe_openmm_top)

        assert openmm_top.getNumAtoms() == gufe_openmm_top.getNumAtoms()
        assert openmm_top.getNumBonds() == gufe_openmm_top.getNumBonds()
        assert openmm_top.getNumChains() == gufe_openmm_top.getNumChains()
        assert openmm_top.getNumResidues() == gufe_openmm_top.getNumResidues()
        assert openmm_top.getNumChains() == gufe_openmm_top.getNumChains()
        # Make sure bond.order is the expected type int or None
        assert all(isinstance(bond.order, (int, type(None))) for bond in openmm_top.bonds())

    @pytest.mark.parametrize("in_pdb_path", ALL_PDB_LOADERS.keys())
    def test_protein_component_openmm_bond(self, in_pdb_path):
        """
        Test that `to_openmm_topology().bonds()` produces a valid bond object.

        We are expecting an StopIteration in order to test that the bonds iterator
        is fully exhausted. Meaning, they are all valid.

        See https://github.com/OpenFreeEnergy/gufe/issues/501
        """
        in_pdb_io = ALL_PDB_LOADERS[in_pdb_path]()
        prot = self.cls.from_pdb_file(in_pdb_io, name="Alice")
        omm_topology_bonds = prot.to_openmm_topology().bonds()

        # Check we are fully exhausting the bond iterator
        with pytest.raises(StopIteration):
            while True:
                repr(next(omm_topology_bonds))

    # Functionality
    def test_eq(self, PDB_181L_path):
        m1 = self.cls.from_pdb_file(PDB_181L_path)
        m2 = self.cls.from_pdb_file(PDB_181L_path)

        assert m1 == m2

    def test_hash_eq(self, PDB_181L_path):
        m1 = self.cls.from_pdb_file(PDB_181L_path)
        m2 = self.cls.from_pdb_file(PDB_181L_path)

        assert hash(m1) == hash(m2)

    def test_neq(self, PDB_181L_path, PDB_181L_mutant):
        m1 = self.cls.from_pdb_file(PDB_181L_path)

        assert m1 != PDB_181L_mutant

    def test_neq_name(self, PDB_181L_path):
        m1 = self.cls.from_pdb_file(PDB_181L_path, name="This")
        m2 = self.cls.from_pdb_file(PDB_181L_path, name="Other")

        assert m1 != m2

    def test_hash_neq(self, PDB_181L_path, PDB_181L_mutant):
        m1 = self.cls.from_pdb_file(PDB_181L_path)

        assert hash(m1) != hash(PDB_181L_mutant)

    def test_hash_neq_name(self, PDB_181L_path):
        m1 = self.cls.from_pdb_file(PDB_181L_path, name="This")
        m2 = self.cls.from_pdb_file(PDB_181L_path, name="Other")

        assert hash(m1) != hash(m2)

    def test_protein_total_charge(self, PDB_181L_path):
        m1 = self.cls.from_pdb_file(PDB_181L_path)

        assert m1.total_charge == 7

    def test_protein_total_charge_thromb(self):
        m1 = self.cls.from_pdb_file(ALL_PDB_LOADERS["thrombin_protein"]())

        assert m1.total_charge == 6

    # these whitespaces in the ion names are intentional to fit pdb formatting
    @pytest.mark.parametrize(
        "ion_name,ion_charge",
        [
            ("SM  ", 3),
            ("Sm  ", 2),
            ("BR  ", -1),
            ("EU3 ", 3),
            (" U4+", 4),
            ("rb  ", 1),
            ("rb2 ", 1),
            ("Hf3 ", 4),
            ("HF  ", 4),
        ],
    )
    def test_pdb_ion_parsing(self, custom_pdb_ion, ion_name, ion_charge):
        pdb_generator = custom_pdb_ion(ion_name)
        pc = self.cls.from_pdb_file(next(pdb_generator))
        expected_total_charge = 8 + ion_charge
        assert pc.total_charge == expected_total_charge

    def test_pdb_ion_invalid(self, custom_pdb_ion):
        ion_name = "ab7 "
        pdb_generator = custom_pdb_ion(ion_name)
        with pytest.raises(ValueError, match="ab7 in residue CL at index 1."):
            self.cls.from_pdb_file(next(pdb_generator))

    def test_protein_valence_error(self):
        """Make sure exception is raised if a valid valence cannot be determined."""
        with pytest.raises(Chem.AtomValenceException, match="Could not set valence of atom id"):
            with mock.patch("rdkit.Chem.rdchem.PeriodicTable.GetValenceList", return_value=[10000000]):
                self.cls.from_pdb_file(ALL_PDB_LOADERS["3tzr_rna"]())


class TestSolvatedPDBComponent(GufeTokenizableTestsMixin, ExplicitMoleculeComponentMixin):
    cls = SolvatedPDBComponent
    repr = "SolvatedPDBComponent(name=Steve)"

    @pytest.fixture(scope="session")
    def instance(self, PDB_a2a_path):
        return self.cls.from_pdb_file(PDB_a2a_path, name="Steve")

    def test_from_pdb_file_sets_box_vectors(self, instance):
        box = instance.box_vectors
        _is_box_shape(box)
        assert _box_vectors_are_in_reduced_form(box)
        assert box[0, 0].m_as(offunit.nanometer) == pytest.approx(6.9587)

    def test_requires_box_vectors(self, PDB_a2a_path):
        prot = ProteinComponent.from_pdb_file(PDB_a2a_path)

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

        comp = factory(path)
        ref = loader(path).topology.getPeriodicBoxVectors()

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
    def test_infer_box_vectors_produces_valid_box(self, factory, path_fixture,
                                                  request):
        path = request.getfixturevalue(path_fixture)

        comp = factory(path, infer_box_vectors=True)
        box = comp.box_vectors

        _is_box_shape(box)
        assert _box_vectors_are_in_reduced_form(box)

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
        comp = self.cls.from_pdbx_file(PDBx_a2a_path, box_vectors=b)
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

        ref = factory(path)
        ref_box = ref.box_vectors

        override = np.eye(3) * 2.0 * offunit.nanometer

        comp = factory(path, box_vectors=override)

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

    def test_cryo_em_dummy_box_raises(self, PDB_181L_path, tmp_path):
        pdb_text = Path(PDB_181L_path).read_text()

        # Insert a CRYST1 record with a 1 Å box
        cryo_em_cryst1 = "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n"
        pdb_path = tmp_path / "cryo_em.pdb"
        pdb_path.write_text(cryo_em_cryst1 + pdb_text)

        pdb = pdbfile.PDBFile(str(pdb_path))

        with pytest.raises(ValueError, match="box_vectors"):
            SolvatedPDBComponent._resolve_box_vectors(pdb)


# class TestProteinMembraneComponent(TestSolvatedPDBComponent):
#     cls = ProteinMembraneComponent
#     repr = "ProteinMembraneComponent(name=Steve)"


def test_no_monomer_info_error(ethane):
    with pytest.raises(TypeError):
        _ = ProteinComponent(rdkit=ethane.to_rdkit())
