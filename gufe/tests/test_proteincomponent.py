# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe
import os
import io
import copy
import importlib.resources
import pytest
from rdkit import Chem
import pooch
import urllib.request
from urllib.error import URLError
from gufe import ProteinComponent

from .test_tokenization import GufeTokenizableTestsMixin

from openmm.app import pdbfile
from openmm import unit
from numpy.testing import assert_almost_equal


try:
    urllib.request.urlopen("https://google.com")
except URLError:
    HAS_INTERNET = False
else:
    HAS_INTERNET = True


_pl_benchmarks = pooch.create(
    path=pooch.os_cache('gufe'),
    base_url="https://github.com/OpenFreeEnergy/openfe-benchmarks/raw/f87f5272cd95d15277de92a7442a9cb71d76615d/openfe_benchmarks/data/",
    version=None,
    registry={
        "cmet_protein.pdb": None,
        "hif2a_protein.pdb": None,
        "mcl1_protein.pdb": None,
        "p38_protein.pdb": None,
        "ptp1b_protein.pdb": None,
        "syk_protein.pdb": None,
        "thrombin_protein.pdb": None,
        "tnsk2_protein.pdb": None,
        "tyk2_protein.pdb": None,
    }
)


@pytest.fixture(
    params=list(_pl_benchmarks.registry.keys()) + ["181l.pdb"]
)
def PLB_PDB_files(request):
    if request.param == '181l.pdb':
        with importlib.resources.path('gufe.tests.data', request.param) as f:
            yield str(f)
    else:
        if not HAS_INTERNET:
            pytest.skip("Offline -- can't fetch PDB file")
        yield _pl_benchmarks.fetch(request.param)


@pytest.fixture
def PDB_181L_mutant(PDB_181L_path):
    # has a single resname flipped to change the resulting sequence
    rdm = Chem.MolFromPDBFile(PDB_181L_path, removeHs=False)

    # important to select the Ca atom - Ca of residue 2
    at = rdm.GetAtomWithIdx(7)
    at.GetMonomerInfo().SetResidueName("X")

    return ProteinComponent.from_rdkit(rdm)


def assert_same_pdb_lines(in_file_path, out_file_path):
    in_lines = []
    if hasattr(in_file_path, "readlines"):
        in_file = in_file_path
    else:
        in_file =  open(in_file_path, "r")

    if isinstance(out_file_path, io.StringIO):
        out_file = out_file_path
        must_close = False
    else:
        out_file = open(out_file_path, mode='r')
        must_close = True

    in_lines = in_file.readlines()
    out_lines = out_file.readlines()

    if must_close:
        out_file.close()

    in_lines = [l for l in in_lines
                if not l.startswith(('REMARK', 'CRYST', '# Created with', 'CONECT'))]
    out_lines = [l for l in out_lines
                 if not l.startswith(('REMARK', 'CRYST', '# Created with', 'CONECT'))]
    in_conect_lines = [l for l in in_lines if l.startswith('CONECT')]
    out_conect_lines = [l for l in in_lines if l.startswith('CONECT')]

    assert in_lines == out_lines


class TestProteinComponent(GufeTokenizableTestsMixin):

    cls = ProteinComponent
    key = "ProteinComponent-72f676dfad8dc1863440a80a53fac542"
    repr = "ProteinComponent(name=Steve)"

    @pytest.fixture
    def instance(self, PDB_181L_OpenMMClean_path):
        return self.cls.from_pdb_file(PDB_181L_OpenMMClean_path, name="Steve")

    def test_from_pdb_file(self, PLB_PDB_files):
        p = self.cls.from_pdb_file(PLB_PDB_files, name="Steve")

        assert isinstance(p, ProteinComponent)
        assert p.name == "Steve"

    def test_from_pdbx_file(self, PDBx_181L_path):
        p = self.cls.from_pdbx_file(str(PDBx_181L_path), name="Steve")

        assert isinstance(p, ProteinComponent)
        assert p.name == "Steve"
        assert p.to_rdkit().GetNumAtoms() == 2639

    def test_from_rdkit(self, PDB_181L_path):
        m = Chem.MolFromPDBFile(PDB_181L_path, removeHs=False)
        p = self.cls.from_rdkit(rdkit=m, name="Steve")

        assert isinstance(p, ProteinComponent)
        assert p.name == "Steve"
        assert p.to_rdkit().GetNumAtoms() == 2639

    # To
    def test_to_rdkit(self, PDB_181L_path):
        pm = self.cls.from_pdb_file(PDB_181L_path)
        rdkitmol = pm.to_rdkit()

        assert isinstance(rdkitmol, Chem.Mol)
        assert rdkitmol.GetNumAtoms() == 2639

    def _test_file_output(self, input_path, output_path, input_type,
                          output_func):
        if input_type == "filename":
            inp = str(output_path)
        elif input_type == "Path":
            inp = output_path
        elif input_type == "StringIO":
            inp = io.StringIO()
        elif input_type == "TextIOWrapper":
            inp = open(output_path, mode='w')

        output_func(inp)

        if input_type == "StringIO":
            inp.seek(0)
            output_path = inp

        if input_type == "TextIOWrapper":
            inp.close()

        assert_same_pdb_lines(in_file_path=str(input_path),
                              out_file_path=output_path)

    @pytest.mark.parametrize('input_type', ['filename', 'Path', 'StringIO',
                                            'TextIOWrapper'])
    def test_to_pdbx_file(self, PDBx_181L_openMMClean_path, tmp_path,
                          input_type):
        p = self.cls.from_pdbx_file(str(PDBx_181L_openMMClean_path), name="Bob")
        out_file_name = "tmp_181L_pdbx.cif"
        out_file = tmp_path / out_file_name

        self._test_file_output(
            input_path=PDBx_181L_openMMClean_path,
            output_path=out_file,
            input_type=input_type,
            output_func=p.to_pdbx_file
        )

    @pytest.mark.parametrize('input_type', ['filename', 'Path', 'StringIO',
                                            'TextIOWrapper'])
    def test_to_pdb_input_types(self, PDB_181L_OpenMMClean_path, tmp_path,
                                input_type):
        p = self.cls.from_pdb_file(str(PDB_181L_OpenMMClean_path), name="Bob")

        self._test_file_output(
            input_path=PDB_181L_OpenMMClean_path,
            output_path=tmp_path / "tmp_181L.pdb",
            input_type=input_type,
            output_func=p.to_pdb_file
        )

    def test_to_pdb_round_trip(self, PLB_PDB_files, tmp_path):
        in_pdb_io = ALL_PDB_LOADERS[in_pdb_path]()

        p = self.cls.from_pdb_file(in_pdb_io, name="Wuff")
        out_file_name = "tmp_"+in_pdb_path+".pdb"
        out_file = tmp_path / out_file_name

        p.to_pdb_file(str(out_file))

        ref_in_pdb_io = ALL_PDB_LOADERS[in_pdb_path]()

        # generate openMM reference file:
        openmm_pdb = pdbfile.PDBFile(ref_in_pdb_io)
        out_ref_file_name = "tmp_"+in_pdb_path+"_openmm_ref.pdb"
        out_ref_file = tmp_path / out_ref_file_name

        pdbfile.PDBFile.writeFile(openmm_pdb.topology, openmm_pdb.positions, file=open(str(out_ref_file), "w"))
        assert_same_pdb_lines(in_file_path=str(out_ref_file), out_file_path=out_file)

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

    def test_to_openmm_positions(self, PLB_PDB_files):
        in_pdb_io = ALL_PDB_LOADERS[in_pdb_path]()
        ref_in_pdb_io =  ALL_PDB_LOADERS[in_pdb_path]()

        openmm_pdb = pdbfile.PDBFile(ref_in_pdb_io)
        openmm_pos = openmm_pdb.positions

        p = self.cls.from_pdb_file(in_pdb_io, name="Bob")
        gufe_openmm_pos = p.to_openmm_positions()

        v1 = gufe_openmm_pos.value_in_unit(unit.nanometer)
        v2 = openmm_pos.value_in_unit(unit.nanometer)

        assert_almost_equal(actual=v1, desired=v2, decimal=6)

    def test_to_openmm_topology(self, PLB_PDB_files):
        in_pdb_io =  ALL_PDB_LOADERS[in_pdb_path]()
        ref_in_pdb_io =  ALL_PDB_LOADERS[in_pdb_path]()

        openmm_pdb = pdbfile.PDBFile(ref_in_pdb_io)
        openmm_top = openmm_pdb.topology

        p = self.cls.from_pdb_file(in_pdb_io, name="Bob")
        gufe_openmm_top = p.to_openmm_topology()
        assert openmm_top.getNumAtoms() == gufe_openmm_top.getNumAtoms()
        for ref_atom, atom in zip(openmm_top.atoms(), gufe_openmm_top.atoms()):
            assert ref_atom.name == atom.name
            assert ref_atom.element == atom.element
            assert ref_atom.index == atom.index
            assert ref_atom.id == atom.id
            assert ref_atom.residue.id == atom.residue.id

        assert openmm_top.getNumResidues() == gufe_openmm_top.getNumResidues()
        for ref_res, res in zip(openmm_top.residues(), gufe_openmm_top.residues()):
            assert len(ref_res) == len(res)
            assert ref_res.name == res.name
            assert ref_res.id == res.id
            assert ref_res.insertionCode == res.insertionCode
            assert ref_res.chain.id == res.chain.id

        for ref_chain, chain in zip(openmm_top.chains(), gufe_openmm_top.chains()):
            assert len(ref_chain) == len(chain)
            assert ref_chain.id == chain.id

        assert openmm_top.getNumBonds() == gufe_openmm_top.getNumBonds()

        # indices are sometimes switched and bonds reordered, but this isn't significant, use set
        ref_bonds = {frozenset({b[0].index, b[1].index}) for b in openmm_top.bonds()}
        other_bonds = {frozenset({b[0].index, b[1].index}) for b in gufe_openmm_top.bonds()}
        assert ref_bonds == other_bonds

        assert openmm_top.getNumAtoms() == gufe_openmm_top.getNumAtoms()
        assert openmm_top.getNumBonds() == gufe_openmm_top.getNumBonds()
        assert openmm_top.getNumChains() == gufe_openmm_top.getNumChains()
        assert openmm_top.getNumResidues() == gufe_openmm_top.getNumResidues()
        assert openmm_top.getNumChains() == gufe_openmm_top.getNumChains()

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
        m1 = self.cls.from_pdb_file(PDB_181L_path, name="This")
        m2 = self.cls.from_pdb_file(PDB_181L_path, name="Other")

        assert m1 != m2

    def test_hash_neq(self, PDB_181L_path, PDB_181L_mutant):
        m1 = self.cls.from_pdb_file(PDB_181L_path)

        assert hash(m1) != hash(PDB_181L_mutant)

    def test_hash_neq_name(self, PDB_181L_path):
        m1 = self.cls.from_pdb_file(PDB_181L_path, name="This")
        m2 = self.cls.from_pdb_file(PDB_181L_path, name="Other")

        assert hash(m1) != hash(m2)

    def test_protein_total_charge(self, PDB_181L_path):
        m1 = self.cls.from_pdb_file(PDB_181L_path)

        assert m1.total_charge == 7

    def test_protein_total_charge_thromb(self):
        m1 = self.cls.from_pdb_file(ALL_PDB_LOADERS["thrombin_protein"]())

        assert m1.total_charge == 6


def test_no_monomer_info_error(ethane):
    with pytest.raises(TypeError):
        _ = ProteinComponent(rdkit=ethane.to_rdkit())
