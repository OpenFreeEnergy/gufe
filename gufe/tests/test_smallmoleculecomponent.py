# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import importlib
import importlib.resources
import openff.toolkit.topology
import os
import hashlib
import pytest

try:
    from openeye import oechem
except ImportError:
    HAS_OECHEM = False
else:
    HAS_OECHEM = oechem.OEChemIsLicensed()
from gufe import SmallMoleculeComponent
from gufe.smallmoleculecomponent import _ensure_ofe_name, _ensure_ofe_version
import gufe
import json
from rdkit import Chem


@pytest.fixture
def alt_ethane():
    return SmallMoleculeComponent(Chem.MolFromSmiles("CC"))


@pytest.fixture
def named_ethane():
    mol = Chem.MolFromSmiles("CC")

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
        assert "SmallMoleculeComponent being renamed" in recwarn[0].message.args[0]
    else:
        assert len(recwarn) == 0

    assert out_name == expected
    assert rdkit.GetProp("ofe-name") == out_name


def test_ensure_ofe_version():
    rdkit = Chem.MolFromSmiles("CC")
    _ensure_ofe_version(rdkit)
    assert rdkit.GetProp("ofe-version") == gufe.__version__


class TestSmallMoleculeComponent:
    def test_rdkit_behavior(self, ethane, alt_ethane):
        # Check that fixture setup is correct (we aren't accidentally
        # testing tautologies)
        assert ethane is not alt_ethane
        assert ethane.to_rdkit() is not alt_ethane.to_rdkit()

    def test_rdkit_independence(self):
        # once we've constructed a Molecule, it is independent from the source
        mol = Chem.MolFromSmiles('CC')
        our_mol = SmallMoleculeComponent.from_rdkit(mol)

        mol.SetProp('foo', 'bar')  # this is the source molecule, not ours
        with pytest.raises(KeyError):
            our_mol.to_rdkit().GetProp('foo')

    def test_rdkit_copy_source_copy(self):
        # we should copy in any properties that were in the source molecule
        mol = Chem.MolFromSmiles('CC')
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

    def test_serialization_cycle(self, named_ethane):
        serialized = named_ethane.to_sdf()
        deserialized = SmallMoleculeComponent.from_sdf_string(serialized)
        reserialized = deserialized.to_sdf()

        assert named_ethane == deserialized
        assert serialized == reserialized

    def test_to_sdf_string(self, named_ethane, serialization_template):
        expected = serialization_template("ethane_template.sdf")
        assert named_ethane.to_sdf() == expected

    def test_from_sdf_string(self, named_ethane, serialization_template):
        sdf_str = serialization_template("ethane_template.sdf")
        assert SmallMoleculeComponent.from_sdf_string(sdf_str) == named_ethane

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
        mol = SmallMoleculeComponent.from_rdkit(rdkit, "ethane")
        assert mol == named_ethane
        assert mol.to_rdkit() is not rdkit

    def test_to_storage_ready(self, named_ethane, serialization_template):
        as_str = serialization_template("ethane_template.sdf")
        bytes_data = as_str.encode("utf-8")
        md5 = hashlib.md5(bytes_data).hexdigest()
        qualname = "gufe.smallmoleculecomponent.SmallMoleculeComponent"
        expected_metadata = {
            ":path:": f"setup/components/{md5}.sdf",
            ":md5:": md5,
            ":class:": qualname,
        }
        expected = {
            named_ethane: gufe.storage.utils.SerializationInfo(
                bytes_data=bytes_data,
                metadata=expected_metadata
            )
        }
        results = named_ethane.to_storage_ready()
        # simple tests so we get errors to tell us what went wrong
        assert len(results) == 1
        assert set(results) == {named_ethane}
        metadata = results[named_ethane].metadata
        assert metadata == expected_metadata
        # full test to ensure that we haven't added fields
        assert results == expected

    def test_from_storage_bytes(self, named_ethane, serialization_template):
        as_str = serialization_template("ethane_template.sdf")
        bytes_data = as_str.encode("utf-8")
        recreated = SmallMoleculeComponent.from_storage_bytes(
            bytes_data,
            load_func=lambda x: None  # no loading required
        )
        assert recreated == named_ethane



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
    sm = SmallMoleculeComponent.from_rdkit(Chem.MolFromSmiles(mol))

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

    def test_bounce_off_file(self, toluene, tmpdir):
        fname = str(tmpdir / 'mol.json')

        with open(fname, 'w') as f:
            f.write(toluene.to_json())
        with open(fname, 'r') as f:
            d = json.load(f)

        assert isinstance(d, dict)
        roundtrip = SmallMoleculeComponent.from_dict(d)

        assert roundtrip == toluene
