# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import pytest
from rdkit import Chem

from gufe import ProteinComponent


@pytest.fixture
def PDB_181L_mutant(PDB_181L_path):
    # has a single resname flipped to change the resulting sequence
    rdm = Chem.MolFromPDBFile(PDB_181L_path, removeHs=False)

    # important to select the Ca atom - Ca of residue 2
    at = rdm.GetAtomWithIdx(7)
    at.GetMonomerInfo().SetResidueName('X')

    return ProteinComponent.from_rdkit(rdm)


def test_from_pdbfile(PDB_181L_path):
    p = ProteinComponent.from_pdbfile(PDB_181L_path, name='Steve')

    assert isinstance(p, ProteinComponent)
    assert p.name == 'Steve'
    assert p._openmm_top.getNumAtoms() == 2639


def test_from_pdbfile_ValueError(PDBx_181L_path):
    with pytest.raises(ValueError):
        _ = ProteinComponent.from_pdbfile(PDBx_181L_path)


@pytest.mark.xfail
def test_from_rdkit(PDB_181L_path):
    m = Chem.MolFromPDBFile(PDB_181L_path, removeHs=False)
    p = ProteinComponent.from_rdkit(m, 'Steve')

    assert isinstance(p, ProteinComponent)
    assert p.name == 'Steve'
    assert p.to_rdkit().GetNumAtoms() == 2639


@pytest.mark.xfail
def test_to_rdkit(PDB_181L_path):
    pm = ProteinComponent.from_pdbfile(PDB_181L_path)
    rdkitmol = pm.to_rdkit()

    assert isinstance(rdkitmol, Chem.Mol)
    assert rdkitmol.GetNumAtoms() == 2639


def test_eq(PDB_181L_path):
    m1 = ProteinComponent.from_pdbfile(PDB_181L_path)
    m2 = ProteinComponent.from_pdbfile(PDB_181L_path)

    assert m1 == m2


def test_hash_eq(PDB_181L_path):
    m1 = ProteinComponent.from_pdbfile(PDB_181L_path)
    m2 = ProteinComponent.from_pdbfile(PDB_181L_path)

    assert hash(m1) == hash(m2)


@pytest.mark.xfail
def test_neq(PDB_181L_path, PDB_181L_mutant):
    m1 = ProteinComponent.from_pdbfile(PDB_181L_path)

    assert m1 != PDB_181L_mutant


def test_neq_name(PDB_181L_path):
    m1 = ProteinComponent.from_pdbfile(PDB_181L_path, name='This')
    m2 = ProteinComponent.from_pdbfile(PDB_181L_path, name='Other')

    assert m1 != m2


@pytest.mark.xfail
def test_hash_neq(PDB_181L_path, PDB_181L_mutant):
    m1 = ProteinComponent.from_pdbfile(PDB_181L_path)

    assert hash(m1) != hash(PDB_181L_mutant)


def test_hash_neq_name(PDB_181L_path):
    m1 = ProteinComponent.from_pdbfile(PDB_181L_path, name='This')
    m2 = ProteinComponent.from_pdbfile(PDB_181L_path, name='Other')

    assert hash(m1) != hash(m2)


@pytest.mark.xfail
def test_protein_total_charge(PDB_181L_path):
    m1 = ProteinComponent.from_pdbfile(PDB_181L_path)

    assert m1.total_charge == 22
