# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import functools
import importlib.resources
import io
import urllib.request
from urllib.error import URLError

import pytest
from openff.units import unit
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


_pl_benchmark_url_pattern = (
    "https://github.com/OpenFreeEnergy/openfe-benchmarks/blob/main/openfe_benchmarks/data/{name}.pdb?raw=true"
)


PDB_BENCHMARK_LOADERS = {
    name: URLFileLike(url=_pl_benchmark_url_pattern.format(name=name)) for name in _benchmark_pdb_names
}

PDB_FILE_LOADERS = {name: lambda: get_test_filename(name) for name in ["181l.pdb"]}

ALL_PDB_LOADERS = dict(**PDB_BENCHMARK_LOADERS, **PDB_FILE_LOADERS)


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
def PDB_181L_OpenMMClean_path():
    with importlib.resources.path("gufe.tests.data", "181l_openmmClean.pdb") as f:
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
            {"solvent": solv_comp, "ligand": ligand}, name=f"{ligand.name}-solvnet"
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
