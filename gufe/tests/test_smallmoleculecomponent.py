# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import importlib
import importlib.resources

try:
    import openff.toolkit.topology
    from openff.units import unit
except ImportError:
    HAS_OFFTK = False
else:
    HAS_OFFTK = True
import json
import logging
import os
from unittest import mock

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

import gufe
from gufe import SmallMoleculeComponent
from gufe.components.explicitmoleculecomponent import _ensure_ofe_name
from gufe.tokenization import TOKENIZABLE_REGISTRY

from .test_explicitmoleculecomponent import ExplicitMoleculeComponentMixin
from .test_tokenization import GufeTokenizableTestsMixin


@pytest.fixture
def alt_ethane():
    mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
    Chem.AllChem.Compute2DCoords(mol)
    return SmallMoleculeComponent(mol)


@pytest.fixture
def named_ethane():
    mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
    Chem.AllChem.Compute2DCoords(mol)
    return SmallMoleculeComponent(mol, name="ethane")


@pytest.mark.parametrize(
    "internal,rdkit_name,name,expected",
    [
        ("", "foo", "", "foo"),
        ("", "", "foo", "foo"),
        ("", "bar", "foo", "foo"),
        ("bar", "", "foo", "foo"),
        ("baz", "bar", "foo", "foo"),
        ("foo", "", "", "foo"),
    ],
)
def test_ensure_ofe_name(internal, rdkit_name, name, expected, recwarn):
    rdkit = Chem.AddHs(Chem.MolFromSmiles("CC"))
    if internal:
        rdkit.SetProp("_Name", internal)

    if rdkit_name:
        rdkit.SetProp("ofe-name", rdkit_name)

    out_name = _ensure_ofe_name(rdkit, name)

    if {rdkit_name, internal} - {"foo", ""}:
        # we should warn if rdkit properties are anything other than 'foo'
        # (expected) or the empty string (not set)
        assert len(recwarn) == 1
        assert "Component being renamed" in recwarn[0].message.args[0]
    else:
        assert len(recwarn) == 0

    assert out_name == expected
    assert rdkit.GetProp("ofe-name") == out_name


class TestSmallMoleculeComponent(GufeTokenizableTestsMixin, ExplicitMoleculeComponentMixin):

    cls = SmallMoleculeComponent
    repr = "SmallMoleculeComponent(name=ethane)"

    @pytest.fixture
    def instance(self, named_ethane):
        return named_ethane

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
        mol = Chem.MolFromSmiles("CC")
        AllChem.Compute2DCoords(mol)
        our_mol = SmallMoleculeComponent.from_rdkit(mol)

        mol.SetProp("foo", "bar")  # this is the source molecule, not ours
        with pytest.raises(KeyError):
            our_mol.to_rdkit().GetProp("foo")

    def test_rdkit_copy_source_copy(self):
        # we should copy in any properties that were in the source molecule
        mol = Chem.MolFromSmiles("CC")
        AllChem.Compute2DCoords(mol)
        mol.SetProp("foo", "bar")
        our_mol = SmallMoleculeComponent.from_rdkit(mol)

        assert our_mol.to_rdkit().GetProp("foo") == "bar"

    def test_equality_and_hash(self, ethane, alt_ethane):
        assert hash(ethane) == hash(alt_ethane)
        assert ethane == alt_ethane

    def test_equality_and_hash_name_differs(self, ethane, named_ethane):
        # names would be used to distinguish different binding modes
        assert hash(ethane) != hash(named_ethane)
        assert ethane != named_ethane

    def test_smiles(self, named_ethane):
        assert named_ethane.smiles == "CC"

    def test_name(self, named_ethane):
        assert named_ethane.name == "ethane"

    def test_empty_name(self, alt_ethane):
        assert alt_ethane.name == ""

    @pytest.mark.xfail
    def test_serialization_cycle(self, named_ethane):
        serialized = named_ethane.to_sdf()
        deserialized = SmallMoleculeComponent.from_sdf_string(serialized)
        reserialized = deserialized.to_sdf()

        assert named_ethane == deserialized
        assert serialized == reserialized

    def test_to_sdf_string(self, named_ethane, ethane_sdf):
        with open(ethane_sdf) as f:
            expected = f.read()

        assert named_ethane.to_sdf() == expected

    @pytest.mark.xfail
    def test_from_sdf_string(self, named_ethane, ethane_sdf):
        with open(ethane_sdf) as f:
            sdf_str = f.read()

        assert SmallMoleculeComponent.from_sdf_string(sdf_str) == named_ethane

    @pytest.mark.xfail
    def test_from_sdf_file(self, named_ethane, ethane_sdf, tmpdir):
        with open(ethane_sdf) as f:
            sdf_str = f.read()
        with open(tmpdir / "temp.sdf", mode="w") as tmpf:
            tmpf.write(sdf_str)

        assert SmallMoleculeComponent.from_sdf_file(tmpdir / "temp.sdf") == named_ethane

    def test_from_sdf_file_junk(self, toluene_mol2_path):
        with pytest.raises(ValueError):
            SmallMoleculeComponent.from_sdf_file(toluene_mol2_path)

    def test_from_sdf_string_multiple_molecules(self, multi_molecule_sdf):
        data = open(multi_molecule_sdf).read()

        with pytest.raises(RuntimeError, match="contains more than 1"):
            SmallMoleculeComponent.from_sdf_string(data)

    def test_from_rdkit(self, named_ethane):
        rdkit = Chem.AddHs(Chem.MolFromSmiles("CC"))
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

    @pytest.mark.parametrize(
        "replace",
        (
            ["name"],
            ["mol"],
            ["name", "mol"],
        ),
    )
    def test_copy_with_replacements(self, named_ethane, replace):
        replacements = {}
        if "name" in replace:
            replacements["name"] = "foo"

        if "mol" in replace:
            # it is a little weird to use copy_with_replacements to replace
            # the whole molecule (possibly keeping the same name), but it
            # should work if someone does! (could more easily imagine only
            # using a new conformer)
            rdmol = Chem.AddHs(Chem.MolFromSmiles("CO"))
            Chem.AllChem.Compute2DCoords(rdmol)
            mol = SmallMoleculeComponent.from_rdkit(rdmol)
            dct = mol._to_dict()
            for item in ["atoms", "bonds", "conformer"]:
                replacements[item] = dct[item]

        new = named_ethane.copy_with_replacements(**replacements)
        if "name" in replace:
            assert new.name == "foo"
        else:
            assert new.name == "ethane"

        if "mol" in replace:
            assert new.smiles == "CO"
        else:
            assert new.smiles == "CC"


@pytest.mark.skipif(not HAS_OFFTK, reason="no openff toolkit available")
class TestSmallMoleculeComponentConversion:
    def test_to_off(self, ethane):
        off_ethane = ethane.to_openff()

        assert isinstance(off_ethane, openff.toolkit.topology.Molecule)

    def test_to_off_name(self, named_ethane):
        off_ethane = named_ethane.to_openff()

        assert off_ethane.name == "ethane"


@pytest.mark.skipif(not HAS_OFFTK, reason="no openff tookit available")
class TestSmallMoleculeComponentPartialCharges:
    @pytest.fixture(scope="function")
    def charged_off_ethane(self, named_ethane):
        off_ethane = named_ethane.to_openff()
        off_ethane.assign_partial_charges(partial_charge_method="am1bcc")
        return off_ethane

    def test_partial_charges_logging(self, charged_off_ethane, caplog):
        caplog.set_level(logging.INFO)
        SmallMoleculeComponent.from_openff(charged_off_ethane)

        assert "Partial charges have been provided" in caplog.text

    def test_partial_charges_zero_warning(self, charged_off_ethane):
        charged_off_ethane.partial_charges[:] = 0 * unit.elementary_charge
        matchmsg = "Partial charges provided all equal to zero"
        with pytest.warns(UserWarning, match=matchmsg):
            SmallMoleculeComponent.from_openff(charged_off_ethane)

    @pytest.mark.parametrize("wrong_charge_val", [-1, 1])
    def test_partial_charges_not_formal_error(self, charged_off_ethane, wrong_charge_val):
        charged_off_ethane.partial_charges[:] = wrong_charge_val * unit.elementary_charge
        with pytest.raises(ValueError, match="Sum of partial charges"):
            SmallMoleculeComponent.from_openff(charged_off_ethane)

    def test_partial_charges_too_few_atoms(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
        Chem.AllChem.Compute2DCoords(mol)
        mol.SetProp("atom.dprop.PartialCharge", "1")

        with pytest.raises(ValueError, match="Incorrect number of"):
            SmallMoleculeComponent.from_rdkit(mol)

    def test_partial_charges_applied_to_atoms(self, caplog):
        """
        Make sure that charges set at the molecule level
        are transferred to atoms and picked up by openFF.
        """
        mol = Chem.AddHs(Chem.MolFromSmiles("C"))
        Chem.AllChem.Compute2DCoords(mol)
        # add some fake charges at the molecule level
        mol.SetProp("atom.dprop.PartialCharge", "-1 0.25 0.25 0.25 0.25")
        caplog.set_level(logging.INFO)

        ofe = SmallMoleculeComponent.from_rdkit(mol)
        assert "Partial charges have been provided" in caplog.text

        # convert to openff and make sure the charges are set
        off_mol = ofe.to_openff()
        assert off_mol.partial_charges is not None
        # check ordering is the same
        rdkit_mol_with_charges = ofe.to_rdkit()
        for i, charge in enumerate(off_mol.partial_charges.m):
            rdkit_atom = rdkit_mol_with_charges.GetAtomWithIdx(i)
            assert rdkit_atom.GetDoubleProp("PartialCharge") == charge

    def test_inconsistent_charges(self, charged_off_ethane):
        """
        An error should be raised if the atom and molecule level
        charges do not match.
        """
        mol = Chem.AddHs(Chem.MolFromSmiles("C"))
        Chem.AllChem.Compute2DCoords(mol)
        # add some fake charges at the molecule level
        mol.SetProp("atom.dprop.PartialCharge", "-1 0.25 0.25 0.25 0.25")
        # set different charges to the atoms
        for atom in mol.GetAtoms():
            atom.SetDoubleProp("PartialCharge", 0)

        # make sure the correct error is raised
        msg = "non-equivalent partial charges between " "atom and molecule properties"
        with pytest.raises(ValueError, match=msg):
            SmallMoleculeComponent.from_rdkit(mol)


@pytest.mark.parametrize(
    "mol, charge",
    [
        ("CC", 0),
        ("CC[O-]", -1),
    ],
)
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

    def test_to_dict_hybridization(self, phenol):
        """
        Make sure dict round trip saves the hybridization
        <https://github.com/OpenFreeEnergy/gufe/issues/407>
        """
        phenol_dict = phenol.to_dict()
        TOKENIZABLE_REGISTRY.clear()
        new_phenol = SmallMoleculeComponent.from_dict(phenol_dict)
        for atom in new_phenol.to_rdkit().GetAtoms():
            if atom.GetAtomicNum() == 6:
                assert atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2

    def test_from_dict_missing_hybridization(self, phenol):
        """
        For backwards compatibility make sure we can create an SMC with missing hybridization info.
        """
        phenol_dict = phenol.to_dict()
        new_atoms = []
        for atom in phenol_dict["atoms"]:
            # remove the hybridization atomic info which should be at index 7
            new_atoms.append(tuple([atom_info for i, atom_info in enumerate(atom) if i != 7]))
        phenol_dict["atoms"] = new_atoms
        with pytest.warns(match="The atom hybridization data was not found and has been set to unspecified."):
            new_phenol = SmallMoleculeComponent.from_dict(phenol_dict)
        # they should be different objects due to the missing hybridization info
        assert new_phenol != phenol
        # make sure the rdkit objects are different
        for atom_hybrid, atom_no_hybrid in zip(phenol.to_rdkit().GetAtoms(), new_phenol.to_rdkit().GetAtoms()):
            assert atom_hybrid.GetHybridization() != atom_no_hybrid.GetHybridization()

    @pytest.mark.skipif(not HAS_OFFTK, reason="no openff toolkit available")
    def test_deserialize_roundtrip(self, toluene, phenol):
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

    @pytest.mark.xfail
    def test_bounce_off_file(self, toluene, tmpdir):
        fname = str(tmpdir / "mol.json")

        with open(fname, "w") as f:
            f.write(toluene.to_json())
        with open(fname) as f:
            d = json.load(f)

        assert isinstance(d, dict)
        roundtrip = SmallMoleculeComponent.from_dict(d)

        assert roundtrip == toluene

    @pytest.mark.skipif(not HAS_OFFTK, reason="no openff toolkit available")
    def test_to_openff_after_serialisation(self, toluene):
        d = toluene.to_dict()

        patch_loc = "gufe.tokenization.TOKENIZABLE_REGISTRY"
        with mock.patch.dict(patch_loc, {}, clear=True):
            t2 = SmallMoleculeComponent.from_dict(d)

        off1 = toluene.to_openff()
        off2 = t2.to_openff()

        assert off1 == off2


@pytest.mark.parametrize("target", ["atom", "bond", "conformer", "mol"])
@pytest.mark.parametrize("dtype", ["int", "bool", "str", "float"])
def test_prop_preservation(ethane, target, dtype):
    # issue 145 make sure props are propagated
    mol = Chem.MolFromSmiles("CC")
    Chem.AllChem.Compute2DCoords(mol)

    if target == "atom":
        obj = mol.GetAtomWithIdx(0)
    elif target == "bond":
        obj = mol.GetBondWithIdx(0)
    elif target == "conformer":
        obj = mol.GetConformer()
    else:
        obj = mol
    if dtype == "int":
        obj.SetIntProp("foo", 1234)
    elif dtype == "bool":
        obj.SetBoolProp("foo", False)
    elif dtype == "str":
        obj.SetProp("foo", "bar")
    elif dtype == "float":
        obj.SetDoubleProp("foo", 1.234)
    else:
        pytest.fail()

    # check that props on rdkit molecules are preserved via to_dict/from_dict cycles
    d = SmallMoleculeComponent(rdkit=mol).to_dict()
    e2 = SmallMoleculeComponent.from_dict(d).to_rdkit()

    if target == "atom":
        obj = e2.GetAtomWithIdx(0)
    elif target == "bond":
        obj = e2.GetBondWithIdx(0)
    elif target == "conformer":
        obj = e2.GetConformer()
    else:
        obj = e2
    if dtype == "int":
        assert obj.GetIntProp("foo") == 1234
    elif dtype == "bool":
        assert obj.GetBoolProp("foo") is False
    elif dtype == "str":
        assert obj.GetProp("foo") == "bar"
    else:
        assert obj.GetDoubleProp("foo") == pytest.approx(1.234)


def test_missing_H_warning():
    m = Chem.MolFromSmiles("CC")
    Chem.AllChem.Compute2DCoords(m)

    with pytest.warns(UserWarning, match="removeHs=False"):
        _ = SmallMoleculeComponent(rdkit=m)


def test_to_openff_5_membered_aromatics():
    # see openfe #278
    m = Chem.MolFromSmiles("c2scnc2C")
    m = Chem.AddHs(m)
    m.Compute2DCoords()

    smc = SmallMoleculeComponent(m)

    offm = smc.to_openff()

    assert offm
