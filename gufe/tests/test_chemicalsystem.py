# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import pytest
import numpy as np

from gufe import ChemicalSystem

from .test_tokenization import GufeTokenizableTestsMixin

def test_ligand_construction(solv_comp, toluene_ligand_comp):
    # sanity checks on construction

    state = ChemicalSystem(
        {'solvent': solv_comp,
         'ligand': toluene_ligand_comp},
    )

    assert len(state.components) == 2
    assert len(state) == 2

    assert list(state) == ['solvent', 'ligand']

    assert state.components['solvent'] == solv_comp
    assert state.components['ligand'] == toluene_ligand_comp
    assert state['solvent'] == solv_comp
    assert state['ligand'] == toluene_ligand_comp


def test_complex_construction(prot_comp, solv_comp, toluene_ligand_comp):
    # sanity checks on construction

    state = ChemicalSystem(
        {'protein': prot_comp,
         'solvent': solv_comp,
         'ligand': toluene_ligand_comp},
    )

    assert len(state.components) == 3
    assert len(state) == 3

    assert list(state) == ['protein', 'solvent', 'ligand']

    assert state.components['protein'] == prot_comp
    assert state.components['solvent'] == solv_comp
    assert state.components['ligand'] == toluene_ligand_comp
    assert state['protein'] == prot_comp
    assert state['solvent'] == solv_comp
    assert state['ligand'] == toluene_ligand_comp


def test_hash_and_eq(prot_comp, solv_comp, toluene_ligand_comp):
    c1 = ChemicalSystem({'protein': prot_comp,
                              'solvent': solv_comp,
                              'ligand': toluene_ligand_comp})

    c2 = ChemicalSystem({'solvent': solv_comp,
                              'ligand': toluene_ligand_comp,
                              'protein': prot_comp})

    assert c1 == c2
    assert hash(c1) == hash(c2)


def test_chemical_system_neq_1(solvated_complex, prot_comp):
    # wrong class
    assert solvated_complex != prot_comp
    assert hash(solvated_complex) != hash(prot_comp)


def test_chemical_system_neq_2(solvated_complex, prot_comp, solv_comp,
                               toluene_ligand_comp):
    # names are different
    complex2 = ChemicalSystem(
        {"protein": prot_comp, "solvent": solv_comp, "ligand": toluene_ligand_comp},
        name="Not quite the same",
    )

    assert solvated_complex != complex2
    assert hash(solvated_complex) != hash(complex2)


def test_chemical_system_neq_3(solvated_complex, prot_comp, solv_comp,
                               toluene_ligand_comp):
    # different unit cell size
    complex2 = ChemicalSystem(
        {'protein': prot_comp,
         'solvent': solv_comp,
         'ligand': toluene_ligand_comp},
        box_vectors=np.array([10, 0, 0] + [np.nan] * 6),
    )
    assert solvated_complex != complex2
    assert hash(solvated_complex) != hash(complex2)


def test_chemical_system_neq_4(solvated_complex, solvated_ligand):
    # different component keys
    assert solvated_complex != solvated_ligand
    assert hash(solvated_complex) != hash(solvated_ligand)


def test_chemical_system_neq_5(solvated_complex, prot_comp, solv_comp,
                               phenol_ligand_comp):
    # same component keys, but different components
    complex2 = ChemicalSystem(
        {'protein': prot_comp,
         'solvent': solv_comp,
         'ligand': phenol_ligand_comp},
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


class TestChemicalSystem(GufeTokenizableTestsMixin):

    cls = ChemicalSystem
    key = "ChemicalSystem-e1cb9ce41e88ee474cf5b962c9388159"

    @pytest.fixture
    def instance(self, solv_comp, toluene_ligand_comp):
        return ChemicalSystem(
            {'solvent': solv_comp,
             'ligand': toluene_ligand_comp},
            box_vectors=np.array([10, 10, 10, 90, 90, 90.]),
            )


@pytest.mark.xfail
class TestChemicalSystemNanBox(GufeTokenizableTestsMixin):

    cls = ChemicalSystem

    @pytest.fixture
    def instance(self, solv_comp, toluene_ligand_comp):
        return ChemicalSystem(
            {'solvent': solv_comp,
             'ligand': toluene_ligand_comp},
            box_vectors=np.array([np.nan, 10, 10, 90, 90, 90.]),
            )
