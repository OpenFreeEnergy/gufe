# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe
import os
import pytest
from rdkit import Chem

from gufe import ProteinComponent

from .test_tokenization import GufeTokenizableTestsMixin

from openmm.app import pdbfile
from openmm import unit
from numpy.testing import assert_almost_equal


@pytest.fixture
def PDB_181L_mutant(PDB_181L_path):
    # has a single resname flipped to change the resulting sequence
    rdm = Chem.MolFromPDBFile(PDB_181L_path, removeHs=False)

    # important to select the Ca atom - Ca of residue 2
    at = rdm.GetAtomWithIdx(7)
    at.GetMonomerInfo().SetResidueName('X')

    return ProteinComponent.from_rdkit(rdm)


def assert_same_pdb_lines(in_file_path, out_file_path):

    in_lines = []
    with  open(in_file_path, "r") as in_pdb:
        in_lines = in_pdb.readlines()
    
    out_lines = []
    with  open(out_file_path, "r") as out_pdb:
        out_lines = out_pdb.readlines()

    
    in_lines = in_lines[1:]
    out_lines = out_lines[1:]
    assert in_lines == out_lines



class TestProteinComponent(GufeTokenizableTestsMixin):

    cls = ProteinComponent

    @pytest.fixture
    def instance(self, PDB_181L_path):
        return self.cls.from_pdb_file(PDB_181L_path, name='Steve')

    # From
    def test_from_pdb_file(self, PDB_181L_path):
        p = self.cls.from_pdb_file(str(PDB_181L_path), name='Steve')

        assert isinstance(p, ProteinComponent)
        assert p.name == 'Steve'
        assert p.to_rdkit().GetNumAtoms() == 2639

    def test_from_pdb_file_benchmark_files(self, PDB_benchmarkFiles):
        for pdb_path in PDB_benchmarkFiles:
            p = self.cls.from_pdb_file(str(pdb_path), name='Steve')

            assert isinstance(p, ProteinComponent)
            assert p.name == 'Steve'

    def test_from_pdb_file_ValueError(self, PDBx_181L_path):
        with pytest.raises(ValueError):
            _ = self.cls.from_pdb_file(PDBx_181L_path)

    
    #@pytest.mark.xfail
    def test_from_pdbx_file(self, PDBx_181L_path):
        p = self.cls.from_pdbx_file(str(PDBx_181L_path), name='Steve')

        assert isinstance(p, ProteinComponent)
        assert p.name == 'Steve'
        assert p.to_rdkit().GetNumAtoms() == 2639

    def test_from_rdkit(self, PDB_181L_path):
        m = Chem.MolFromPDBFile(PDB_181L_path, removeHs=False)
        p = self.cls.from_rdkit(rdkit=m, name='Steve')

        assert isinstance(p, ProteinComponent)
        assert p.name == 'Steve'
        assert p.to_rdkit().GetNumAtoms() == 2639

    # To
    def test_to_rdkit(self, PDB_181L_path):
        pm = self.cls.from_pdb_file(PDB_181L_path)
        rdkitmol = pm.to_rdkit()

        assert isinstance(rdkitmol, Chem.Mol)
        assert rdkitmol.GetNumAtoms() == 2639

    def test_to_pdbx_file(self, PDB_181L_path, tmp_path):
        p = self.cls.from_pdb_file(str(PDB_181L_path), name='Bob')
        out_file_name = "tmp_181L_pdbx.cif"
        out_file = tmp_path / out_file_name

        p.to_pdbx_file(str(out_file))        
        with  open(str(out_file), "w") as out_file2:
            p.to_pdbx_file(out_file2)
        

    def test_to_pdb_file(self, PDB_181L_path, tmp_path):
        p = self.cls.from_pdb_file(str(PDB_181L_path), name='Wuff')
        out_file_name = "tmp_181L_pdb.pdb"
        out_file = tmp_path / out_file_name

        p.to_pdb_file(str(out_file))
        
        with  open(str(out_file), "w") as out_file2:           
            p.to_pdb_file(out_file2)

    def test_io_pdb_comparison(self, PDB_181L_OpenMMClean_path, tmp_path):
        out_file_name = "tmp_" + os.path.basename(PDB_181L_OpenMMClean_path)
        out_file = tmp_path / out_file_name

        p = self.cls.from_pdb_file(PDB_181L_OpenMMClean_path, name="Bob")
        _ = p.to_pdb_file(str(out_file))

        assert_same_pdb_lines(PDB_181L_OpenMMClean_path, str(out_file))

    def test_dummy_from_dict(self, PDB_181L_OpenMMClean_path):
        p = self.cls.from_pdb_file(PDB_181L_OpenMMClean_path, name="Bob")
        gufe_dict = p.to_dict()
        p2 = self.cls.from_dict(gufe_dict)

        assert p == p2

    def test_to_openmm_positions(self, PDB_181L_OpenMMClean_path):
        openmm_pdb = pdbfile.PDBFile(open(PDB_181L_OpenMMClean_path, "r"))
        openmm_pos = openmm_pdb.positions

        p = self.cls.from_pdb_file(PDB_181L_OpenMMClean_path, name="Bob")
        gufe_openmm_pos = p.to_openmm_positions()

        v1 = gufe_openmm_pos.value_in_unit(unit.nanometer)
        v2 = openmm_pos.value_in_unit(unit.nanometer)

        assert_almost_equal(actual=v1, desired=v2, decimal=6)

    def test_to_openmm_positions_bench(self, PDB_benchmarkFiles):
        for pdb_path in PDB_benchmarkFiles:
            openmm_pdb = pdbfile.PDBFile(open(pdb_path, "r"))
            openmm_pos = openmm_pdb.positions

            p = self.cls.from_pdb_file(pdb_path, name="BobThePudel")
            gufe_openmm_pos = p.to_openmm_positions()

            v1 = gufe_openmm_pos.value_in_unit(unit.nanometer)
            v2 = openmm_pos.value_in_unit(unit.nanometer)

            assert_almost_equal(actual=v1, desired=v2, decimal=6)

    def test_to_openmm_topology(self, PDB_181L_OpenMMClean_path, tmp_path):
        openmm_pdb = pdbfile.PDBFile(open(PDB_181L_OpenMMClean_path, "r"))
        openmm_top = openmm_pdb.topology

        p = self.cls.from_pdb_file(PDB_181L_OpenMMClean_path, name="Bob")
        gufe_openmm_top = p.to_openmm_topology()

        assert openmm_top.getNumAtoms() == gufe_openmm_top.getNumAtoms()
        assert openmm_top.getNumBonds() == gufe_openmm_top.getNumBonds()
        assert openmm_top.getNumChains() == gufe_openmm_top.getNumChains()
        assert openmm_top.getNumResidues() == gufe_openmm_top.getNumResidues()
        assert openmm_top.getNumChains() == gufe_openmm_top.getNumChains()

        assert gufe_openmm_top.getPeriodicBoxVectors() is None
        assert gufe_openmm_top.getUnitCellDimensions() is None

    def test_to_openmm_topology_bench(self, PDB_benchmarkFiles):
        for in_pdb_path in PDB_benchmarkFiles:

            openmm_pdb = pdbfile.PDBFile(open(in_pdb_path, "r"))
            openmm_top = openmm_pdb.topology

            p = self.cls.from_pdb_file(in_pdb_path, name="Bob")
            gufe_openmm_top = p.to_openmm_topology()

            assert openmm_top.getNumAtoms() == gufe_openmm_top.getNumAtoms()
            assert openmm_top.getNumBonds() == gufe_openmm_top.getNumBonds()
            assert openmm_top.getNumChains() == gufe_openmm_top.getNumChains()
            assert openmm_top.getNumChains() == gufe_openmm_top.getNumChains()

            nresis1 = openmm_top.getNumResidues()
            nresi2 = gufe_openmm_top.getNumResidues()
            assert nresis1 == nresi2

            pbvs1 = gufe_openmm_top.getPeriodicBoxVectors()
            pbvs2 = openmm_top.getPeriodicBoxVectors()
            v1 = pbvs1.value_in_unit(unit.nanometer)
            v2 = pbvs2.value_in_unit(unit.nanometer)
            assert_almost_equal(
                actual=v1,
                desired=v2,
                decimal=6,
                err_msg="the pbcVs are not equal")

            uD = gufe_openmm_top.getUnitCellDimensions()
            v1 = uD.value_in_unit(unit.nanometer)

            uD2 = openmm_top.getUnitCellDimensions()
            v2 = uD2.value_in_unit(unit.nanometer)

            assert_almost_equal(
                actual=v1,
                desired=v2,
                decimal=6,
                err_msg="the unitcellDims are not equal")

    def test_io_pdb_comparison_bench(self, PDB_benchmarkFiles, tmp_path):
        failures = []
        for in_pdb_path in PDB_benchmarkFiles:
            out_file_name = "tmp_" + os.path.basename(in_pdb_path)
            out_file = tmp_path / out_file_name

            try:
                p = self.cls.from_pdb_file(in_pdb_path, name="bench")
                _ = p.to_pdb_file(str(out_file))

                assert_same_pdb_lines(in_pdb_path, str(out_file))
            except KeyError as err:
                failures.append((in_pdb_path, err))

        if(len(failures) > 0):
            raise Exception("IFailed: " + str(failures))
    # Functionality

    def test_eq(self, PDB_181L_path):
        m1 = self.cls.from_pdb_file(PDB_181L_path)
        m2 = self.cls.from_pdb_file(PDB_181L_path)

        assert m1 == m2

    def test_hash_eq(self, PDB_181L_path):
        m1 = self.cls.from_pdb_file(PDB_181L_path)
        m2 = self.cls.from_pdb_file(PDB_181L_path)

        assert hash(m1) == hash(m2)

    def test_neq(self, PDB_181L_path, PDB_181L_mutant):
        m1 = self.cls.from_pdb_file(PDB_181L_path)

        assert m1 != PDB_181L_mutant

    def test_neq_name(self, PDB_181L_path):
        m1 = self.cls.from_pdb_file(PDB_181L_path, name='This')
        m2 = self.cls.from_pdb_file(PDB_181L_path, name='Other')

        assert m1 != m2

    def test_hash_neq(self, PDB_181L_path, PDB_181L_mutant):
        m1 = self.cls.from_pdb_file(PDB_181L_path)

        assert hash(m1) != hash(PDB_181L_mutant)

    def test_hash_neq_name(self, PDB_181L_path):
        m1 = self.cls.from_pdb_file(PDB_181L_path, name='This')
        m2 = self.cls.from_pdb_file(PDB_181L_path, name='Other')

        assert hash(m1) != hash(m2)

    def test_protein_total_charge(self, PDB_181L_path):
        m1 = self.cls.from_pdb_file(PDB_181L_path)

        assert m1.total_charge == 7

    def test_protein_total_charge_thromb(self, PDB_thrombin_path):
        m1 = self.cls.from_pdb_file(PDB_thrombin_path)

        assert m1.total_charge == 6
