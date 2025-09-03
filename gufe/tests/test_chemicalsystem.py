# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import numpy as np
import pytest

from gufe import ChemicalSystem, SmallMoleculeComponent, SolventComponent
from gufe.components import ProteinComponent

from .test_tokenization import GufeTokenizableTestsMixin
from ..components.explicitmoleculecomponent import ExplicitMoleculeComponent


def test_ligand_construction(solv_comp, toluene_ligand_comp):
    # sanity checks on construction

    state = ChemicalSystem(
        {"solvent": solv_comp, "ligand": toluene_ligand_comp},
    )

    assert len(state.components) == 2
    assert len(state) == 2

    assert list(state) == ["solvent", "ligand"]

    assert state.components["solvent"] == solv_comp
    assert state.components["ligand"] == toluene_ligand_comp
    assert state["solvent"] == solv_comp
    assert state["ligand"] == toluene_ligand_comp


def test_complex_construction(prot_comp, solv_comp, toluene_ligand_comp):
    # sanity checks on construction

    state = ChemicalSystem(
        {"protein": prot_comp, "solvent": solv_comp, "ligand": toluene_ligand_comp},
    )

    assert len(state.components) == 3
    assert len(state) == 3

    assert list(state) == ["protein", "solvent", "ligand"]

    assert state.components["protein"] == prot_comp
    assert state.components["solvent"] == solv_comp
    assert state.components["ligand"] == toluene_ligand_comp
    assert state["protein"] == prot_comp
    assert state["solvent"] == solv_comp
    assert state["ligand"] == toluene_ligand_comp


def test_hash_and_eq(prot_comp, solv_comp, toluene_ligand_comp):
    c1 = ChemicalSystem({"protein": prot_comp, "solvent": solv_comp, "ligand": toluene_ligand_comp})

    c2 = ChemicalSystem({"solvent": solv_comp, "ligand": toluene_ligand_comp, "protein": prot_comp})

    assert c1 == c2
    assert hash(c1) == hash(c2)


def test_chemical_system_component_diff(solvated_complex, solvated_ligand):
    comps_diff = solvated_complex.component_diff(solvated_ligand)

    assert isinstance(comps_diff, tuple)

    comps_complex, comps_ligand = comps_diff

    assert isinstance(comps_complex, tuple)
    assert isinstance(comps_ligand, tuple)

    assert len(comps_complex) == 1
    assert len(comps_ligand) == 0

    assert isinstance(comps_complex[0], ProteinComponent)

    # A.component_diff(B) equals the reversed output of B.component_diff(A)
    assert solvated_complex.component_diff(solvated_ligand) == solvated_ligand.component_diff(solvated_complex)[::-1]


def test_chemical_system_component_diff_incompatible_comparison(solvated_complex, phenol_ligand_comp):
    with pytest.raises(
        TypeError, match=r"`other` must be an instance of `ChemicalSystem`, not `SmallMoleculeComponent`"
    ):
        solvated_complex.component_diff(phenol_ligand_comp)


def test_chemical_system_neq_1(solvated_complex, prot_comp):
    # wrong class
    assert solvated_complex != prot_comp
    assert hash(solvated_complex) != hash(prot_comp)


def test_chemical_system_neq_2(solvated_complex, prot_comp, solv_comp, toluene_ligand_comp):
    # names are different
    complex2 = ChemicalSystem(
        {"protein": prot_comp, "solvent": solv_comp, "ligand": toluene_ligand_comp},
        name="Not quite the same",
    )

    assert solvated_complex != complex2
    assert hash(solvated_complex) != hash(complex2)


def test_chemical_system_neq_4(solvated_complex, solvated_ligand):
    # different component keys
    assert solvated_complex != solvated_ligand
    assert hash(solvated_complex) != hash(solvated_ligand)


def test_chemical_system_neq_5(solvated_complex, prot_comp, solv_comp, phenol_ligand_comp):
    # same component keys, but different components
    complex2 = ChemicalSystem(
        {"protein": prot_comp, "solvent": solv_comp, "ligand": phenol_ligand_comp},
    )
    assert solvated_complex != complex2
    assert hash(solvated_complex) != hash(complex2)


@pytest.mark.xfail
def test_complex_system_charge(solvated_complex):
    # protein = 22, ligand = 0, solvent = 0
    assert solvated_complex.total_charge == 22


def test_ligand_system_charge(solvated_ligand):
    assert solvated_ligand.total_charge == 0


def test_sorting(solvated_complex, solvated_ligand):
    order1 = [solvated_complex, solvated_ligand, solvated_ligand]
    order2 = [solvated_ligand, solvated_complex, solvated_ligand]

    assert sorted(order1) == sorted(order2)


def test_isin_wrong_type(solvated_complex):
    with pytest.raises(TypeError, match="`item` must be an instance or subclass of `Component`"):
        solvated_complex.isin(float)

def test_isin_instance(solvated_complex, prot_comp, toluene_ligand_comp, phenol_ligand_comp):
    # check for present instances don't return matches
    assert solvated_complex.isin(prot_comp) is True

def test_isin_type(solvated_complex):
    # check for present types don't return matches
    assert solvated_complex.isin(ProteinComponent) is True

def test_isin_instance_return_matches(solvated_complex, prot_comp, toluene_ligand_comp, phenol_ligand_comp):
    # check for present instances
    isin, matches = solvated_complex.isin(prot_comp, return_matches=True)
    assert isin is True
    assert matches == [prot_comp]

    isin, matches = solvated_complex.isin(toluene_ligand_comp, return_matches=True)
    assert isin is True
    assert matches == [toluene_ligand_comp]

    # check for absent instance
    isin, matches = solvated_complex.isin(phenol_ligand_comp, return_matches=True)
    assert isin is False
    assert matches == []

def test_isin_type_return_matches(solvated_ligand):
    # check for present types
    isin, matches = solvated_ligand.isin(ProteinComponent, return_matches=True)
    assert isin is False
    assert matches == []

    isin, matches = solvated_ligand.isin(SmallMoleculeComponent, return_matches=True)
    assert isin is True
    assert matches == [solvated_ligand.components["ligand"]]

    isin, matches = solvated_ligand.isin(SolventComponent, return_matches=True)
    assert isin is True
    assert matches == [solvated_ligand.components["solvent"]]



class TestChemicalSystem(GufeTokenizableTestsMixin):
    cls = ChemicalSystem
    repr = "ChemicalSystem(name=, components={'solvent': SolventComponent(name=O, K+, Cl-), 'ligand': SmallMoleculeComponent(name=toluene)})"

    @pytest.fixture
    def instance(self, solv_comp, toluene_ligand_comp):
        return ChemicalSystem(
            {"solvent": solv_comp, "ligand": toluene_ligand_comp},
        )
