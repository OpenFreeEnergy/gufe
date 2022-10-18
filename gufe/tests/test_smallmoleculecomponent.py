# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import importlib
import importlib.resources
import openff.toolkit.topology
import os
import pytest

try:
    from openeye import oechem
except ImportError:
    HAS_OECHEM = False
else:
    HAS_OECHEM = oechem.OEChemIsLicensed()
from gufe import SmallMoleculeComponent
from gufe.components.explicitmoleculecomponent import (
    _ensure_ofe_name, _ensure_ofe_version
)
import gufe
import json
from rdkit import Chem
from rdkit.Chem import AllChem
from gufe.tokenization import TOKENIZABLE_REGISTRY

from .test_tokenization import GufeTokenizableTestsMixin

@pytest.fixture
def alt_ethane():
    mol = Chem.MolFromSmiles("CC")
    Chem.AllChem.Compute2DCoords(mol)
    return SmallMoleculeComponent(mol)

@pytest.fixture
def named_ethane():
    mol = Chem.MolFromSmiles("CC")
    Chem.AllChem.Compute2DCoords(mol)
    return SmallMoleculeComponent(mol, name='ethane')


@pytest.mark.parametrize('internal,rdkit_name,name,expected', [
    ('', 'foo', '', 'foo'),
    ('', '', 'foo', 'foo'),
    ('', 'bar', 'foo', 'foo'),
    ('bar', '', 'foo', 'foo'),
    ('baz', 'bar', 'foo', 'foo'),
    ('foo', '', '', 'foo'),
])
def test_ensure_ofe_name(internal, rdkit_name, name, expected, recwarn):
    rdkit = Chem.MolFromSmiles("CC")
    if internal:
        rdkit.SetProp('_Name', internal)

    if rdkit_name:
        rdkit.SetProp('ofe-name', rdkit_name)

    out_name = _ensure_ofe_name(rdkit, name)

    if {rdkit_name, internal} - {'foo', ''}:
        # we should warn if rdkit properties are anything other than 'foo'
        # (expected) or the empty string (not set)
        assert len(recwarn) == 1
        assert "Component being renamed" in recwarn[0].message.args[0]
    else:
        assert len(recwarn) == 0

    assert out_name == expected
    assert rdkit.GetProp("ofe-name") == out_name


def test_ensure_ofe_version():
    rdkit = Chem.MolFromSmiles("CC")
    _ensure_ofe_version(rdkit)
    assert rdkit.GetProp("ofe-version") == gufe.__version__


class TestSmallMoleculeComponent(GufeTokenizableTestsMixin):

    cls = SmallMoleculeComponent
    key = "SmallMoleculeComponent-3a1b343b46ec93300bc74d83c133637a"

    @pytest.fixture
    def instance(self, named_ethane):
        return named_ethane

    def test_rdkit_behavior(self, ethane, alt_ethane):
        # Check that fixture setup is correct (we aren't accidentally
        # testing tautologies)
        assert ethane is not alt_ethane
        assert ethane.to_rdkit() is not alt_ethane.to_rdkit()


    def test_error_missing_conformers(self):
        mol = Chem.MolFromSmiles("CC")
        with pytest.raises(ValueError, match="conformer"):
            SmallMoleculeComponent(mol)

    def test_warn_multiple_conformers(self):
        mol = Chem.MolFromSmiles("CC")
        AllChem.EmbedMultipleConfs(mol)
        with pytest.warns(UserWarning, match="conformers. Only"):
            SmallMoleculeComponent(mol)

    def test_rdkit_independence(self):
        # once we've constructed a Molecule, it is independent from the source
        mol = Chem.MolFromSmiles('CC')
        AllChem.Compute2DCoords(mol)
        our_mol = SmallMoleculeComponent.from_rdkit(mol)

        mol.SetProp('foo', 'bar')  # this is the source molecule, not ours
        with pytest.raises(KeyError):
            our_mol.to_rdkit().GetProp('foo')

    def test_rdkit_copy_source_copy(self):
        # we should copy in any properties that were in the source molecule
        mol = Chem.MolFromSmiles('CC')
        AllChem.Compute2DCoords(mol)
        mol.SetProp('foo', 'bar')
        our_mol = SmallMoleculeComponent.from_rdkit(mol)

        assert our_mol.to_rdkit().GetProp('foo') == 'bar'

    def test_equality_and_hash(self, ethane, alt_ethane):
        assert hash(ethane) == hash(alt_ethane)
        assert ethane == alt_ethane

    def test_equality_and_hash_name_differs(self, ethane, named_ethane):
        # names would be used to distinguish different binding modes
        assert hash(ethane) != hash(named_ethane)
        assert ethane != named_ethane

    def test_smiles(self, named_ethane):
        assert named_ethane.smiles == 'CC'

    def test_name(self, named_ethane):
        assert named_ethane.name == 'ethane'

    def test_empty_name(self, alt_ethane):
        assert alt_ethane.name == ''

    @pytest.mark.xfail
    def test_serialization_cycle(self, named_ethane):
        serialized = named_ethane.to_sdf()
        deserialized = SmallMoleculeComponent.from_sdf_string(serialized)
        reserialized = deserialized.to_sdf()

        assert named_ethane == deserialized
        assert serialized == reserialized

    def test_to_sdf_string(self, named_ethane, serialization_template):
        expected = serialization_template("ethane_template.sdf")
        assert named_ethane.to_sdf() == expected

    @pytest.mark.xfail
    def test_from_sdf_string(self, named_ethane, serialization_template):
        sdf_str = serialization_template("ethane_template.sdf")
        assert SmallMoleculeComponent.from_sdf_string(sdf_str) == named_ethane

    @pytest.mark.xfail
    def test_from_sdf_file(self, named_ethane, serialization_template,
                           tmpdir):
        sdf_str = serialization_template("ethane_template.sdf")
        with open(tmpdir / "temp.sdf", mode='w') as tmpf:
            tmpf.write(sdf_str)

        assert SmallMoleculeComponent.from_sdf_file(tmpdir / "temp.sdf") == named_ethane

    def test_from_sdf_file_junk(self, toluene_mol2_path):
        with pytest.raises(ValueError):
            SmallMoleculeComponent.from_sdf_file(toluene_mol2_path)

    def test_from_sdf_string_multiple_molecules(self, multi_molecule_sdf):
        data = open(multi_molecule_sdf, 'r').read()

        with pytest.raises(RuntimeError, match="contains more than 1"):
            SmallMoleculeComponent.from_sdf_string(data)

    def test_from_rdkit(self, named_ethane):
        rdkit = Chem.MolFromSmiles("CC")
        AllChem.Compute2DCoords(rdkit)
        mol = SmallMoleculeComponent.from_rdkit(rdkit, "ethane")
        assert mol == named_ethane
        assert mol.to_rdkit() is not rdkit

    def test_serialization_cycle_smiles(self, named_ethane):
        # check a regression against the smiles changing on serialization
        dct = named_ethane.to_dict()
        TOKENIZABLE_REGISTRY.clear()
        copy = SmallMoleculeComponent.from_dict(dct)
        assert named_ethane == copy
        assert named_ethane is not copy
        assert named_ethane.smiles == copy.smiles


class TestSmallMoleculeComponentConversion:
    def test_to_off(self, ethane):
        off_ethane = ethane.to_openff()

        assert isinstance(off_ethane, openff.toolkit.topology.Molecule)

    def test_to_off_name(self, named_ethane):
        off_ethane = named_ethane.to_openff()

        assert off_ethane.name == 'ethane'

    @pytest.mark.skipif(not HAS_OECHEM, reason="No OEChem available")
    def test_to_oechem(self, ethane):
        oec_ethane = ethane.to_openeye()

        assert isinstance(oec_ethane, oechem.OEMol)


@pytest.mark.parametrize('mol, charge', [
    ('CC', 0), ('CC[O-]', -1),
])
def test_total_charge_neutral(mol, charge):
    mol = Chem.MolFromSmiles(mol)
    AllChem.Compute2DCoords(mol)
    sm = SmallMoleculeComponent.from_rdkit(mol)

    assert sm.total_charge == charge


def test_sorting(ethane, alt_ethane):
    # to check larger containers it's important to be able to sort
    order1 = [ethane, alt_ethane, ethane]
    order2 = [alt_ethane, ethane, ethane]

    assert sorted(order1) == sorted(order2)


class TestSmallMoleculeSerialization:
    def test_to_dict(self, phenol):
        d = phenol.to_dict()

        assert isinstance(d, dict)

    def test_deserialize_roundtrip(self, toluene, phenol):
        # TODO: Currently roundtripping via openff adds Hydrogens even when
        #       they weren't in the original input molecule.
        roundtrip = SmallMoleculeComponent.from_dict(phenol.to_dict())

        assert roundtrip == phenol
        assert roundtrip != toluene

        # check the coordinates, these aren't in the hash
        # but are important to preserve
        pos1 = phenol.to_openff().conformers
        pos2 = roundtrip.to_openff().conformers
        assert len(pos1) == len(pos2)
        for x, y in zip(pos1, pos2):
            assert (x == y).all()

    # TODO: determine if we want to add our own serializers for e.g. JSON
    # based on `to_dict`
    @pytest.mark.xfail
    def test_bounce_off_file(self, toluene, tmpdir):
        fname = str(tmpdir / 'mol.json')

        with open(fname, 'w') as f:
            f.write(toluene.to_json())
        with open(fname, 'r') as f:
            d = json.load(f)

        assert isinstance(d, dict)
        roundtrip = SmallMoleculeComponent.from_dict(d)

        assert roundtrip == toluene
