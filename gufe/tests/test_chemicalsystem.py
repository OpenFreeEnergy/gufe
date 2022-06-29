# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import pytest
import numpy as np

import json
import hashlib

import gufe


def test_ligand_construction(solv_comp, toluene_ligand_comp):
    # sanity checks on construction

    state = gufe.ChemicalSystem(
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

    state = gufe.ChemicalSystem(
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
    c1 = gufe.ChemicalSystem({'protein': prot_comp,
                              'solvent': solv_comp,
                              'ligand': toluene_ligand_comp})

    c2 = gufe.ChemicalSystem({'solvent': solv_comp,
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
    complex2 = gufe.ChemicalSystem(
        {"protein": prot_comp, "solvent": solv_comp, "ligand": toluene_ligand_comp},
        name="Not quite the same",
    )

    assert solvated_complex != complex2
    assert hash(solvated_complex) != hash(complex2)


def test_chemical_system_neq_3(solvated_complex, prot_comp, solv_comp,
                               toluene_ligand_comp):
    # different unit cell size
    complex2 = gufe.ChemicalSystem(
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
    complex2 = gufe.ChemicalSystem(
        {'protein': prot_comp,
         'solvent': solv_comp,
         'ligand': phenol_ligand_comp},
    )
    assert solvated_complex != complex2
    assert hash(solvated_complex) != hash(complex2)


def test_chemical_system_to_storage_ready(solvated_ligand):
    storage_ready = solvated_ligand.to_storage_ready()
    solvent = solvated_ligand['solvent']
    ligand = solvated_ligand['ligand']
    assert len(storage_ready) == 3
    assert set(storage_ready.keys()) == {solvated_ligand, solvent, ligand}

    # TODO: it would be really nice to provide convenience methods for
    # testing some of the stuff before. A lot of this is transferable to
    # other problems.
    info = storage_ready[solvated_ligand]
    dct = json.loads(info.bytes_data.decode('utf-8'))
    assert set(dct) == {'components', 'box_vectors', 'name'}

    assert dct['name'] is None

    for val in dct['box_vectors']:
        assert np.isnan(val)

    dct_ligand = dct['components']['ligand']
    dct_solvent = dct['components']['solvent']
    # ensure that we've wrapped these correctly
    assert set(dct_ligand) == {":path:", ":class:", ":module:", ":md5:"}
    assert set(dct_solvent) == {":path:", ":class:", ":module:", ":md5:"}

    # NOTE: this is weird. I'm still having trouble getting arrays of NaN to
    # be recognized as default values.
    expected_non_defaults = ['components']
    expected_hash_dict = {key: dct[key] for key in expected_non_defaults}
    expected_hash = hashlib.md5(
        json.dumps(expected_hash_dict, sort_keys=True).encode('utf-8')
    ).hexdigest()
    assert info.metadata[":md5:"] == expected_hash


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
