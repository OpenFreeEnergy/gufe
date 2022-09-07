# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe
import os
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


def line_by_line_comparison(in_file_path, out_file_path):
    
        in_pdb = open(in_file_path, "r")
        out_pdb = open(out_file_path, "r")
        
        in_lines = in_pdb.readlines()
        out_lines = out_pdb.readlines()
        
        not_equal = []
        assert len(in_lines) == len(out_lines)

        for i, in_line in enumerate(in_lines):
            out_line = out_lines[i]
            if(in_line == out_line):
                continue
            else:
                not_equal.append((in_line, out_line))          
            
        not_equal_found = len(not_equal)>0
        if(not_equal_found): 
            print("not Equals") 
            print("\n".join(map(lambda x: "in: "+x[0]+"\nout: "+x[1], not_equal))) 
            
        
        if(not_equal_found):
            return False    
        else: 
            return True               

class TestProteinComponent(GufeTokenizableTestsMixin):

    cls = ProteinComponent

    @pytest.fixture
    def instance(self, PDB_181L_path):
        return self.cls.from_pdbfile(PDB_181L_path, name='Steve')

    # From 
    def test_from_pdbfile(self, PDB_181L_path):
        p = self.cls.from_pdbfile(str(PDB_181L_path), name='Steve')

        assert isinstance(p, ProteinComponent)
        assert p.name == 'Steve'
        assert p.to_rdkit().GetNumAtoms() == 2639

    def test_from_pdbfile_benchmark_files(self, PDB_benchmarkFiles):
        for pdb_path in PDB_benchmarkFiles:
            p = self.cls.from_pdbfile(str(pdb_path), name='Steve')

            assert isinstance(p, ProteinComponent)
            assert p.name == 'Steve'
       
    def test_from_pdbfile_ValueError(self, PDBx_181L_path):
        with pytest.raises(ValueError):
            _ = self.cls.from_pdbfile(PDBx_181L_path)

    @pytest.mark.xfail
    def test_from_pdbxfile(self, PDBx_181L_path):
        p = self.cls.from_pdbxfile(str(PDBx_181L_path), name='Steve')

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
        pm = self.cls.from_pdbfile(PDB_181L_path)
        rdkitmol = pm.to_rdkit()

        assert isinstance(rdkitmol, Chem.Mol)
        assert rdkitmol.GetNumAtoms() == 2639

    def test_to_pdbxfile(self, PDB_181L_path, tmp_path):
        p = self.cls.from_pdbfile(str(PDB_181L_path), name='Bob')
        out_path_prefix = "tmp_181L_pdbx.cif"
        out_file = tmp_path / out_path_prefix
        p.to_pdbxFile(str(out_file))
    
    def test_to_pdbfile(self, PDB_181L_path, tmp_path):
        p = self.cls.from_pdbfile(str(PDB_181L_path), name='Wuff')
        out_path_prefix = "tmp_181L_pdb.pdb"
        out_file = tmp_path / out_path_prefix
        p.to_pdbFile(str(out_file))
    
    def test_io_pdb_comparison(self, PDB_181L_OpenMMClean_path, tmp_path):
        out_path_prefix = "tmp_"+os.path.basename(PDB_181L_OpenMMClean_path)
        out_file = tmp_path / out_path_prefix
        
        p = self.cls.from_pdbfile(PDB_181L_OpenMMClean_path, name="Bob")
        _ = p.to_pdbFile(str(out_file))            
    
        assert line_by_line_comparison(PDB_181L_OpenMMClean_path, str(out_file))    
    
    def test_to_openmm(self, PDB_181L_OpenMMClean_path, tmp_path):
        from openmm.app import pdbfile
        openmm_pdb = pdbfile.PDBFile(open(PDB_181L_OpenMMClean_path, "r"))
        openmm_top = openmm_pdb.topology
        
        p = self.cls.from_pdbfile(PDB_181L_OpenMMClean_path, name="Bob")
        gufe_openmm_top = p.to_openmm_topology()
        
        assert openmm_top.getNumAtoms() == gufe_openmm_top.getNumAtoms()
        assert openmm_top.getNumBonds() == gufe_openmm_top.getNumBonds()
        assert openmm_top.getNumChains() == gufe_openmm_top.getNumChains()
        assert openmm_top.getNumResidues() == gufe_openmm_top.getNumResidues()
        assert openmm_top.getNumChains() == gufe_openmm_top.getNumChains()
        assert openmm_top.getPeriodicBoxVectors() == gufe_openmm_top.getPeriodicBoxVectors() 
        assert openmm_top.getUnitCellDimensions() == gufe_openmm_top.getUnitCellDimensions() 


    def test_to_openmm_bench(self, PDB_benchmarkFiles):
        from openmm.app import pdbfile
        for in_pdb_path in PDB_benchmarkFiles:

            openmm_pdb = pdbfile.PDBFile(open(in_pdb_path, "r"))
            openmm_top = openmm_pdb.topology
            
            p = self.cls.from_pdbfile(in_pdb_path, name="Bob")
            gufe_openmm_top = p.to_openmm_topology()
            
            # assert openmm_top == gufe_openmm_top
            assert openmm_top.getNumAtoms() == gufe_openmm_top.getNumAtoms()
            assert openmm_top.getNumBonds() == gufe_openmm_top.getNumBonds()
            assert openmm_top.getNumChains() == gufe_openmm_top.getNumChains()
            assert openmm_top.getNumResidues() == gufe_openmm_top.getNumResidues()
            assert openmm_top.getNumChains() == gufe_openmm_top.getNumChains()
            assert openmm_top.getPeriodicBoxVectors() == gufe_openmm_top.getPeriodicBoxVectors() 
            assert openmm_top.getUnitCellDimensions() == gufe_openmm_top.getUnitCellDimensions() 

        
    def test_io_pdb_comparison_bench(self, PDB_benchmarkFiles,tmp_path):
        failures = []
        for in_pdb_path in PDB_benchmarkFiles:
            out_path_prefix = "tmp_"+os.path.basename(in_pdb_path)
            out_file = tmp_path / out_path_prefix
            
            try:
                p = self.cls.from_pdbfile(in_pdb_path, name="bench")
                _ = p.to_pdbFile(str(out_file))            
                
                assert line_by_line_comparison(in_pdb_path, str(out_file))  
            except KeyError as err:
                failures.append((in_pdb_path, err)) 
        
        if(len(failures)>0):
            raise Exception("IFailed: "+str(failures))

    # Functionality
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
        m1 = self.cls.from_pdbfile(PDB_181L_path)

        assert m1.total_charge == 7
        
    def test_protein_total_charge_thromb(self, PDB_thrombin_path):
        m1 = self.cls.from_pdbfile(PDB_thrombin_path)

        assert m1.total_charge == 6
