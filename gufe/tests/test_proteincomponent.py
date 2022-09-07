# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import pytest
from rdkit import Chem

from gufe import ProteinComponent

from .test_tokenize import GufeTokenizableTestsMixin


@pytest.fixture
def PDB_181L_mutant(PDB_181L_path):
    # has a single resname flipped to change the resulting sequence
    rdm = Chem.MolFromPDBFile(PDB_181L_path, removeHs=False)

    # important to select the Ca atom - Ca of residue 2
    at = rdm.GetAtomWithIdx(7)
    at.GetMonomerInfo().SetResidueName('X')
    
    return ProteinComponent.from_rdkit(rdm)




class TestProteinComponent(GufeTokenizableTestsMixin):

    cls = ProteinComponent

    @pytest.fixture
    def instance(self, PDB_181L_path):
        return self.cls.from_pdbfile(PDB_181L_path, name='Steve')


    def test_from_pdbfile(self, PDB_181L_path):
        print(str(PDB_181L_path))
        p = self.cls.from_pdbfile(str(PDB_181L_path), name='Steve')

        assert isinstance(p, ProteinComponent)
        assert p.name == 'Steve'
        assert p.to_rdkit().GetNumAtoms() == 2639


    def test_from_pdbfile_ValueError(self, PDBx_181L_path):
        with pytest.raises(ValueError):
            _ = self.cls.from_pdbfile(PDBx_181L_path)


    def test_from_rdkit(self, PDB_181L_path):
        m = Chem.MolFromPDBFile(PDB_181L_path, removeHs=False)
        p = self.cls.from_rdkit(rdkit=m, name='Steve')

        assert isinstance(p, ProteinComponent)
        assert p.name == 'Steve'
        assert p.to_rdkit().GetNumAtoms() == 2639


    def test_to_rdkit(self, PDB_181L_path):
        pm = self.cls.from_pdbfile(PDB_181L_path)
        rdkitmol = pm.to_rdkit()

        assert isinstance(rdkitmol, Chem.Mol)
        assert rdkitmol.GetNumAtoms() == 2639


    def test_eq(self, PDB_181L_path):
        m1 = self.cls.from_pdbfile(PDB_181L_path)
        m2 = self.cls.from_pdbfile(PDB_181L_path)

        assert m1 == m2


    def test_hash_eq(self, PDB_181L_path):
        m1 = self.cls.from_pdbfile(PDB_181L_path)
        m2 = self.cls.from_pdbfile(PDB_181L_path)

        assert hash(m1) == hash(m2)


    def test_neq(self, PDB_181L_path, PDB_181L_mutant):
        m1 = self.cls.from_pdbfile(PDB_181L_path)

        assert m1 != PDB_181L_mutant


    def test_neq_name(self, PDB_181L_path):
        m1 = self.cls.from_pdbfile(PDB_181L_path, name='This')
        m2 = self.cls.from_pdbfile(PDB_181L_path, name='Other')

        assert m1 != m2


    def test_hash_neq(self, PDB_181L_path, PDB_181L_mutant):
        m1 = self.cls.from_pdbfile(PDB_181L_path)

        assert hash(m1) != hash(PDB_181L_mutant)


    def test_hash_neq_name(self, PDB_181L_path):
        m1 = self.cls.from_pdbfile(PDB_181L_path, name='This')
        m2 = self.cls.from_pdbfile(PDB_181L_path, name='Other')

        assert hash(m1) != hash(m2)


    def test_protein_total_charge(self, PDB_181L_path):
        print("IN: ", str(PDB_181L_path))
        m1 = self.cls.from_pdbfile(PDB_181L_path)

        assert m1.total_charge == 7 #Todo: verify, used to be 22
        
        
    def test_protein_total_charge_thromb(self, PDB_thrombin_path):
        print("IN: ", str(PDB_thrombin_path))
        m1 = self.cls.from_pdbfile(PDB_thrombin_path)

        assert m1.total_charge == -2
