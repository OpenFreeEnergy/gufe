# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import importlib
import json
import pathlib

import numpy as np
import pytest
from openff.units import unit
from rdkit import Chem

import gufe
from gufe import LigandAtomMapping, SmallMoleculeComponent

from .test_tokenization import GufeTokenizableTestsMixin


def mol_from_smiles(smiles: str) -> gufe.SmallMoleculeComponent:
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    m.Compute2DCoords()

    return gufe.SmallMoleculeComponent(m)


@pytest.fixture(scope="session")
def simple_mapping():
    """Disappearing oxygen on end

    C C O

    C C
    """
    molA = mol_from_smiles("CCO")
    molB = mol_from_smiles("CC")

    m = LigandAtomMapping(molA, molB, componentA_to_componentB={0: 0, 1: 1})

    return m


@pytest.fixture(scope="session")
def other_mapping():
    """Disappearing middle carbon

    C C O

    C   C
    """
    molA = mol_from_smiles("CCO")
    molB = mol_from_smiles("CC")

    m = LigandAtomMapping(molA, molB, componentA_to_componentB={0: 0, 2: 1})

    return m


@pytest.fixture
def annotated_simple_mapping(simple_mapping):
    mapping = LigandAtomMapping(
        simple_mapping.componentA,
        simple_mapping.componentB,
        simple_mapping.componentA_to_componentB,
        annotations={"foo": "bar"},
    )
    return mapping


@pytest.fixture(scope="session")
def benzene_maps():
    MAPS = {
        "phenol": {
            0: 0,
            1: 1,
            2: 2,
            3: 3,
            4: 4,
            5: 5,
            6: 6,
            7: 7,
            8: 8,
            9: 9,
            10: 12,
            11: 11,
        },
        "anisole": {
            0: 5,
            1: 6,
            2: 7,
            3: 8,
            4: 9,
            5: 10,
            6: 11,
            7: 12,
            8: 13,
            9: 14,
            10: 2,
            11: 15,
        },
    }
    return MAPS


@pytest.fixture(scope="session")
def benzene_phenol_mapping(benzene_transforms, benzene_maps):
    molA = SmallMoleculeComponent(benzene_transforms["benzene"].to_rdkit())
    molB = SmallMoleculeComponent(benzene_transforms["phenol"].to_rdkit())
    m = LigandAtomMapping(molA, molB, benzene_maps["phenol"])
    return m


@pytest.fixture(scope="session")
def benzene_anisole_mapping(benzene_transforms, benzene_maps):
    molA = SmallMoleculeComponent(benzene_transforms["benzene"].to_rdkit())
    molB = SmallMoleculeComponent(benzene_transforms["anisole"].to_rdkit())
    m = LigandAtomMapping(molA, molB, benzene_maps["anisole"])
    return m


@pytest.fixture(scope="session")
def atom_mapping_basic_test_files():
    # a dict of {filenames.strip(mol2): SmallMoleculeComponent} for a simple
    # set of ligands
    files = {}
    for f in [
        "1,3,7-trimethylnaphthalene",
        "1-butyl-4-methylbenzene",
        "2,6-dimethylnaphthalene",
        "2-methyl-6-propylnaphthalene",
        "2-methylnaphthalene",
        "2-naftanol",
        "methylcyclohexane",
        "toluene",
    ]:
        with importlib.resources.path("gufe.tests.data.lomap_basic", f + ".mol2") as fn:
            mol = Chem.MolFromMol2File(str(fn), removeHs=False)
            files[f] = SmallMoleculeComponent(mol, name=f)

    return files


def test_atommapping_usage(simple_mapping):
    assert simple_mapping.componentA_to_componentB[1] == 1
    assert simple_mapping.componentA_to_componentB.get(2, None) is None
    assert simple_mapping.annotations == {}

    with pytest.raises(KeyError):
        simple_mapping.componentA_to_componentB[3]


def test_mapping_inversion(benzene_phenol_mapping):
    assert benzene_phenol_mapping.componentB_to_componentA == {
        0: 0,
        1: 1,
        2: 2,
        3: 3,
        4: 4,
        5: 5,
        6: 6,
        7: 7,
        8: 8,
        9: 9,
        11: 11,
        12: 10,
    }


def test_mapping_distances(benzene_phenol_mapping):
    d = benzene_phenol_mapping.get_distances()

    ref = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.34005502, 0.0]

    assert isinstance(d, np.ndarray)
    for i, r in zip(d, ref):
        assert i == pytest.approx(r)


def test_uniques(atom_mapping_basic_test_files):
    mapping = LigandAtomMapping(
        componentA=atom_mapping_basic_test_files["methylcyclohexane"],
        componentB=atom_mapping_basic_test_files["toluene"],
        componentA_to_componentB={0: 6, 1: 7, 2: 8, 3: 9, 4: 10, 5: 11, 6: 12},
    )

    assert list(mapping.componentA_unique) == [
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
    ]
    assert list(mapping.componentB_unique) == [0, 1, 2, 3, 4, 5, 13, 14]


def test_modification(benzene_phenol_mapping):
    # check that we get a copy of the mapping and we can't modify
    AtoB = benzene_phenol_mapping.componentA_to_componentB
    before = len(AtoB)

    AtoB.pop(10)

    assert len(benzene_phenol_mapping.componentA_to_componentB) == before


def test_atommapping_hash(simple_mapping, other_mapping):
    # these two mappings map the same molecules, but with a different mapping
    assert simple_mapping is not other_mapping


def test_draw_mapping_cairo(tmpdir, simple_mapping):
    with tmpdir.as_cwd():
        simple_mapping.draw_to_file("test.png")
        filed = pathlib.Path("test.png")
        assert filed.exists()


def test_draw_mapping_svg(tmpdir, other_mapping):
    with tmpdir.as_cwd():
        d2d = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(600, 300, 300, 300)
        other_mapping.draw_to_file("test.svg", d2d=d2d)
        filed = pathlib.Path("test.svg")
        assert filed.exists()


class TestLigandAtomMappingSerialization:
    def test_deserialize_roundtrip(self, benzene_phenol_mapping, benzene_anisole_mapping):
        roundtrip = LigandAtomMapping.from_dict(benzene_phenol_mapping.to_dict())

        assert roundtrip == benzene_phenol_mapping

        # We don't check coordinates since that's already done in guefe for
        # SmallMoleculeComponent

        assert roundtrip != benzene_anisole_mapping

    def test_file_roundtrip(self, benzene_phenol_mapping, tmpdir):
        with tmpdir.as_cwd():
            with open("tmpfile.json", "w") as f:
                f.write(json.dumps(benzene_phenol_mapping.to_dict()))

            with open("tmpfile.json") as f:
                d = json.load(f)

            assert isinstance(d, dict)
            roundtrip = LigandAtomMapping.from_dict(d)

            assert roundtrip == benzene_phenol_mapping


def test_annotated_atommapping_hash_eq(simple_mapping, annotated_simple_mapping):
    assert annotated_simple_mapping != simple_mapping
    assert hash(annotated_simple_mapping) != hash(simple_mapping)


def test_annotation_immutability(annotated_simple_mapping):
    annot1 = annotated_simple_mapping.annotations
    annot1["foo"] = "baz"
    annot2 = annotated_simple_mapping.annotations
    assert annot1 != annot2
    assert annot2 == {"foo": "bar"}


def test_with_annotations(simple_mapping, annotated_simple_mapping):
    new_annot = simple_mapping.with_annotations({"foo": "bar"})
    assert new_annot == annotated_simple_mapping


def test_with_fancy_annotations(simple_mapping):
    m = simple_mapping.with_annotations({"thing": 4.0 * unit.nanometer})

    assert m.key

    m2 = LigandAtomMapping.from_dict(m.to_dict())

    assert m == m2


class TestLigandAtomMappingBoundsChecks:
    @pytest.fixture
    def molA(self):
        # 9 atoms
        return mol_from_smiles("CCO")

    @pytest.fixture
    def molB(self):
        # 11 atoms
        return mol_from_smiles("CCC")

    def test_too_large_A(self, molA, molB):
        with pytest.raises(ValueError, match="invalid index for ComponentA"):
            LigandAtomMapping(componentA=molA, componentB=molB, componentA_to_componentB={9: 5})

    def test_too_small_A(self, molA, molB):
        with pytest.raises(ValueError, match="invalid index for ComponentA"):
            LigandAtomMapping(componentA=molA, componentB=molB, componentA_to_componentB={-2: 5})

    def test_too_large_B(self, molA, molB):
        with pytest.raises(ValueError, match="invalid index for ComponentB"):
            LigandAtomMapping(componentA=molA, componentB=molB, componentA_to_componentB={5: 11})

    def test_too_small_B(self, molA, molB):
        with pytest.raises(ValueError, match="invalid index for ComponentB"):
            LigandAtomMapping(componentA=molA, componentB=molB, componentA_to_componentB={5: -1})


class TestLigandAtomMapping(GufeTokenizableTestsMixin):
    cls = LigandAtomMapping
    repr = "LigandAtomMapping(componentA=SmallMoleculeComponent(name=), componentB=SmallMoleculeComponent(name=), componentA_to_componentB={0: 0, 1: 1}, annotations={'foo': 'bar'})"

    @pytest.fixture
    def instance(self, annotated_simple_mapping):
        return annotated_simple_mapping

    def test_id_key(self, instance):
        i2 = self.cls.from_dict(instance.to_dict())

        assert instance.key == i2.key

    def test_keyed_dict(self, instance):
        i2 = self.cls.from_dict(instance.to_dict())

        assert instance.to_keyed_dict() == i2.to_keyed_dict()
