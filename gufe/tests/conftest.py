# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import importlib.resources
import urllib.request
from urllib.error import URLError
import io
import functools
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

import gufe

try:
    urllib.request.urlopen("https://google.com")
except URLError:
    HAS_INTERNET = False
else:
    HAS_INTERNET = True


## helper functions


def load_url_data(url, encoding='utf-8'):
    if not HAS_INTERNET:
        pytest.skip("Skipping because internet seems faulty")

    req = urllib.request.urlopen(url)
    # we convert to StringIO because req.read() returns bytes, not string
    return io.StringIO(req.read().decode(encoding))


def get_test_filename(filename):
    with importlib.resources.path('gufe.tests.data', filename) as file:
        return str(file)


## PDBs for input/output testing


_benchmark_pdb_names = [
        "cmet_protein",
        "hif2a_protein",
        "mcl1_protein",
        "p38_protein",
        "ptp1b_protein",
        "syk_protein",
        "thrombin_protein",
        "tnsk2_protein",
        "tyk2_protein",
        ]


_pl_benchmark_url_pattern = (
    "https://github.com/OpenFreeEnergy/openfe-benchmarks/blob/main/openfe_benchmarks/data/{name}.pdb?raw=true"
)


PDB_BENCHMARK_LOADERS = {
    name: functools.partial(
        load_url_data,
        url=_pl_benchmark_url_pattern.format(name=name)
    )
    for name in _benchmark_pdb_names
}

PDB_FILE_LOADERS = {
    name: lambda: get_test_filename(name)
    for name in ["181l.pdb"]
}

ALL_PDB_LOADERS = dict(**PDB_BENCHMARK_LOADERS, **PDB_FILE_LOADERS)


## data file paths

@pytest.fixture
def serialization_template():
    def inner(filename):
        loc = "gufe.tests.data"
        tmpl = importlib.resources.read_text(loc, filename)
        return tmpl.format(OFE_VERSION=gufe.__version__)

    return inner


@pytest.fixture
def toluene_mol2_path():
    with importlib.resources.path('gufe.tests.data', 'toluene.mol2') as f:
        yield str(f)


@pytest.fixture
def multi_molecule_sdf():
    fn = 'multi_molecule.sdf'
    with importlib.resources.path('gufe.tests.data', fn) as f:
        yield str(f)


@pytest.fixture
def PDB_181L_path():
    with importlib.resources.path('gufe.tests.data', '181l.pdb') as f:
        yield str(f)


@pytest.fixture
def PDB_181L_OpenMMClean_path():
    with importlib.resources.path('gufe.tests.data',
                                  '181l_openmmClean.pdb') as f:
        yield str(f)


@pytest.fixture
def PDB_thrombin_path():
    with importlib.resources.path('gufe.tests.data',
                                  'thrombin_protein.pdb') as f:
        yield str(f)


@pytest.fixture
def PDBx_181L_path():
    with importlib.resources.path('gufe.tests.data',
                                  '181l.cif') as f:
        yield str(f)

@pytest.fixture
def PDBx_181L_openMMClean_path():
    with importlib.resources.path('gufe.tests.data',
                                  '181l_openmmClean.cif') as f:
        yield str(f)

## RDKit molecules

@pytest.fixture
def benzene_modifications():
    with importlib.resources.path('gufe.tests.data',
                                  'benzene_modifications.sdf') as f:
        supp = Chem.SDMolSupplier(str(f), removeHs=False)

        mols = list(supp)

    return {m.GetProp('_Name'): m for m in mols}


## Components


@pytest.fixture
def benzene(benzene_modifications):
    return gufe.SmallMoleculeComponent(benzene_modifications["benzene"])


@pytest.fixture
def toluene(benzene_modifications):
    return gufe.SmallMoleculeComponent(benzene_modifications["toluene"])


@pytest.fixture
def phenol(benzene_modifications):
    return gufe.SmallMoleculeComponent(benzene_modifications['phenol'])


@pytest.fixture
def benzonitrile(benzene_modifications):
    return gufe.SmallMoleculeComponent(benzene_modifications["benzonitrile"])


@pytest.fixture
def anisole(benzene_modifications):
    return gufe.SmallMoleculeComponent(benzene_modifications["anisole"])


@pytest.fixture
def benzaldehyde(benzene_modifications):
    return gufe.SmallMoleculeComponent(benzene_modifications["benzaldehyde"])


@pytest.fixture
def styrene(benzene_modifications):
    return gufe.SmallMoleculeComponent(benzene_modifications["styrene"])


@pytest.fixture(scope="session")
def ethane():
    mol = Chem.MolFromSmiles("CC")
    AllChem.Compute2DCoords(mol)
    return gufe.SmallMoleculeComponent(mol)


@pytest.fixture
def prot_comp(PDB_181L_path):
    yield gufe.ProteinComponent.from_pdb_file(PDB_181L_path)


@pytest.fixture
def solv_comp():
    yield gufe.SolventComponent(positive_ion="K", negative_ion="Cl")


@pytest.fixture
def toluene_ligand_comp(benzene_modifications):
    yield gufe.SmallMoleculeComponent.from_rdkit(benzene_modifications["toluene"])


@pytest.fixture
def phenol_ligand_comp(benzene_modifications):
    yield gufe.SmallMoleculeComponent.from_rdkit(benzene_modifications["phenol"])


## ChemicalSystems


@pytest.fixture
def solvated_complex(prot_comp, solv_comp, toluene_ligand_comp):
    return gufe.ChemicalSystem(
        {"protein": prot_comp, "solvent": solv_comp, "ligand": toluene_ligand_comp},
    )


@pytest.fixture
def solvated_ligand(solv_comp, toluene_ligand_comp):
    return gufe.ChemicalSystem(
        {"ligand": toluene_ligand_comp, "solvent": solv_comp},
    )

@pytest.fixture
def vacuum_ligand(toluene_ligand_comp):
    return gufe.ChemicalSystem(
        {"ligand": toluene_ligand_comp},
    )
