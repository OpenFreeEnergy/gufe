# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import pytest
from openff.toolkit.topology import Molecule

import gufe


@pytest.fixture
def prot_comp(PDB_181L_path):
    yield gufe.ProteinComponent.from_pdbfile(PDB_181L_path)


@pytest.fixture
def solv_comp():
    yield gufe.SolventComponent(ions=('K', 'Cl'))


@pytest.fixture
def lig_comp(benzene_modifications):
    yield gufe.SmallMoleculeComponent.from_rdkit(benzene_modifications['toluene'])


def test_construction(prot_comp, solv_comp, lig_comp):
    # sanity checks on construction

    state = gufe.ChemicalState(
        {'protein': prot_comp,
         'solvent': solv_comp,
         'ligand': lig_comp},
    )

    assert len(state.components) == 3

    assert state.components['protein'] == prot_comp
    assert state.components['solvent'] == solv_comp
    assert state.components['ligand'] == lig_comp


def test_hash_and_eq(prot_comp, solv_comp, lig_comp):
    c1 = gufe.ChemicalState({'protein': prot_comp,
                             'solvent': solv_comp,
                             'ligand': lig_comp})

    c2 = gufe.ChemicalState({'solvent': solv_comp,
                             'ligand': lig_comp,
                             'protein': prot_comp})

    assert c1 == c2
    assert hash(c1) == hash(c2)
