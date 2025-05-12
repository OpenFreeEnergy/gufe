# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import functools
import importlib.resources
import io
import urllib.request
from urllib.error import URLError

import pooch
import pytest
from openff.units import unit
from packaging.version import Version
from rdkit import Chem
from rdkit.Chem import AllChem

import gufe
from gufe.tests.test_protocol import DummyProtocol

try:
    urllib.request.urlopen("https://google.com")
except URLError:
    HAS_INTERNET = False
else:
    HAS_INTERNET = True

try:
    from openmm import Platform

    OPENMM_VERSION = Version(Platform.getOpenMMVersion())
# Support OpenMM as a soft dep
except ModuleNotFoundError:
    OPENMM_VERSION = None


class URLFileLike:
    def __init__(self, url, encoding="utf-8"):
        self.url = url
        self.encoding = encoding
        self.data = None

    def __call__(self):
        if not HAS_INTERNET:  # pragma: no-cover
            pytest.skip("Skipping because internet seems faulty")

        if self.data is None:
            req = urllib.request.urlopen(self.url)
            self.data = req.read().decode(self.encoding)

        return io.StringIO(self.data)


def get_test_filename(filename):
    with importlib.resources.path("gufe.tests.data", filename) as file:
        return str(file)


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


_pl_benchmark_url_pattern = "https://github.com/OpenFreeEnergy/openfe-benchmarks/tree/f577a88d94b6deae1c0de7eb926edf48dd42b72d/openfe_benchmarks/data/{name}.pdb?raw=true"


# This will duplcate the PDB_BENCHMARK_LOADERS
# We can keep PDB_BENCHMARK_LOADERS arround just in case we want to use it again someday
POOCH_CACHE = pooch.os_cache("gufe")
PDB_FILE = pooch.create(
    path=POOCH_CACHE,
    base_url="doi:10.5281/zenodo.15376306",
    registry={
        "hif2a_protein.pdb": "md5:158711011c6b85b1d55ad68a559ca07b",
        "cmet_protein.pdb": "md5:0a38968fe9c09a49da2f1c57e84196f0",
        "mcl1_protein.pdb": "md5:a2e57c14a925ee0dd9f158eb7ad5f413",
        "p38_protein.pdb": "md5:52e1ba80c73c7f710219580bddb71de3",
        "ptp1b_protein.pdb": "md5:f754c3aa3ea1617d99a4dd21fb158967",
        "syk_protein.pdb": "md5:0c631f4b88ed7bbd35ee660214426cd2",
        "thrombin_protein.pdb": "md5:261b8f040b389188e8c0cf14fbab5775",
        "tnsk2_protein.pdb": "md5:aa13bef540d061ed66ef4240886a39ec",
        "tyk2_protein.pdb": "md5:48d447290ee637ce8e3255cfa572297a",
    },
    retry_if_failed=10,
)


# TODO We can either fixtureize this or make it lazy to speed things up
class PoochFileLike:
    """Wrapper so we can use the calls that exist to PDB_BENCHMARK_LOADERS"""

    def __init__(self, name):
        self.name = name

    def __call__(self):
        return PDB_FILE.fetch(f"{self.name}.pdb")


PDB_ZENODO_LOADERS = {name: PoochFileLike(name) for name in _benchmark_pdb_names}


# PDB_BENCHMARK_LOADERS = {
#    name: URLFileLike(url=_pl_benchmark_url_pattern.format(name=name)) for name in _benchmark_pdb_names
# }

PDB_FILE_LOADERS = {name: lambda: get_test_filename(name) for name in ["181l.pdb"]}

ALL_PDB_LOADERS = dict(**PDB_ZENODO_LOADERS, **PDB_FILE_LOADERS)


@pytest.fixture
def ethane_sdf():
    with importlib.resources.path("gufe.tests.data", "ethane.sdf") as f:
        yield str(f)


@pytest.fixture
def toluene_mol2_path():
    with importlib.resources.path("gufe.tests.data", "toluene.mol2") as f:
        yield str(f)


@pytest.fixture
def multi_molecule_sdf():
    fn = "multi_molecule.sdf"
    with importlib.resources.path("gufe.tests.data", fn) as f:
        yield str(f)


@pytest.fixture
def PDB_181L_path():
    with importlib.resources.path("gufe.tests.data", "181l.pdb") as f:
        yield str(f)


@pytest.fixture
def offxml_settings_path():
    with importlib.resources.path("gufe.tests.data", "offxml_settings.json") as f:
        yield str(f)


@pytest.fixture
def all_settings_path():
    with importlib.resources.path("gufe.tests.data", "all_settings.json") as f:
        yield str(f)


@pytest.fixture
def PDB_thrombin_path():
    with importlib.resources.path("gufe.tests.data", "thrombin_protein.pdb") as f:
        yield str(f)


@pytest.fixture
def PDBx_181L_path():
    with importlib.resources.path("gufe.tests.data", "181l.cif") as f:
        yield str(f)


@pytest.fixture
def PDBx_181L_openMMClean_path():
    with importlib.resources.path("gufe.tests.data", "181l_openmmClean.cif") as f:
        yield str(f)


@pytest.fixture(scope="session")
def benzene_modifications():
    with importlib.resources.path("gufe.tests.data", "benzene_modifications.sdf") as f:
        supp = Chem.SDMolSupplier(str(f), removeHs=False)

        mols = list(supp)

    return {m.GetProp("_Name"): m for m in mols}


@pytest.fixture(scope="session")
def benzene_transforms(benzene_modifications):
    return {k: gufe.SmallMoleculeComponent(v) for k, v in benzene_modifications.items()}


@pytest.fixture
def benzene(benzene_modifications):
    return gufe.SmallMoleculeComponent(benzene_modifications["benzene"])


@pytest.fixture
def toluene(benzene_modifications):
    return gufe.SmallMoleculeComponent(benzene_modifications["toluene"])


@pytest.fixture
def phenol(benzene_modifications):
    return gufe.SmallMoleculeComponent(benzene_modifications["phenol"])


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
    mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
    AllChem.Compute2DCoords(mol)
    return gufe.SmallMoleculeComponent(mol)


@pytest.fixture
def prot_comp(PDB_181L_path):
    yield gufe.ProteinComponent.from_pdb_file(PDB_181L_path)


@pytest.fixture
def solv_comp():
    yield gufe.SolventComponent(positive_ion="K", negative_ion="Cl", ion_concentration=0.0 * unit.molar)


@pytest.fixture
def toluene_ligand_comp(benzene_modifications):
    yield gufe.SmallMoleculeComponent.from_rdkit(benzene_modifications["toluene"])


@pytest.fixture
def phenol_ligand_comp(benzene_modifications):
    yield gufe.SmallMoleculeComponent.from_rdkit(benzene_modifications["phenol"])


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


@pytest.fixture
def absolute_transformation(solvated_ligand, solvated_complex):
    return gufe.Transformation(
        solvated_ligand,
        solvated_complex,
        protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
        mapping=None,
    )


@pytest.fixture
def complex_equilibrium(solvated_complex):
    return gufe.NonTransformation(
        solvated_complex,
        protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
    )


@pytest.fixture
def benzene_variants_star_map_transformations(
    benzene,
    toluene,
    phenol,
    benzonitrile,
    anisole,
    benzaldehyde,
    styrene,
    prot_comp,
    solv_comp,
):

    variants = [toluene, phenol, benzonitrile, anisole, benzaldehyde, styrene]

    # define the solvent chemical systems and transformations between
    # benzene and the others
    solvated_ligands = {}
    solvated_ligand_transformations = {}

    solvated_ligands["benzene"] = gufe.ChemicalSystem({"solvent": solv_comp, "ligand": benzene}, name="benzene-solvent")

    for ligand in variants:
        solvated_ligands[ligand.name] = gufe.ChemicalSystem(
            {"solvent": solv_comp, "ligand": ligand}, name=f"{ligand.name}-solvent"
        )
        solvated_ligand_transformations[("benzene", ligand.name)] = gufe.Transformation(
            solvated_ligands["benzene"],
            solvated_ligands[ligand.name],
            protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
            mapping=None,
        )

    # define the complex chemical systems and transformations between
    # benzene and the others
    solvated_complexes = {}
    solvated_complex_transformations = {}

    solvated_complexes["benzene"] = gufe.ChemicalSystem(
        {"protein": prot_comp, "solvent": solv_comp, "ligand": benzene},
        name="benzene-complex",
    )

    for ligand in variants:
        solvated_complexes[ligand.name] = gufe.ChemicalSystem(
            {"protein": prot_comp, "solvent": solv_comp, "ligand": ligand},
            name=f"{ligand.name}-complex",
        )
        solvated_complex_transformations[("benzene", ligand.name)] = gufe.Transformation(
            solvated_complexes["benzene"],
            solvated_complexes[ligand.name],
            protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
            mapping=None,
        )

    return list(solvated_ligand_transformations.values()), list(solvated_complex_transformations.values())


@pytest.fixture
def benzene_variants_star_map(benzene_variants_star_map_transformations):
    solvated_ligand_transformations, solvated_complex_transformations = benzene_variants_star_map_transformations
    return gufe.AlchemicalNetwork(solvated_ligand_transformations + solvated_complex_transformations)


@pytest.fixture
def benzene_variants_ligand_star_map(benzene_variants_star_map_transformations):
    solvated_ligand_transformations, _ = benzene_variants_star_map_transformations
    return gufe.AlchemicalNetwork(solvated_ligand_transformations)
