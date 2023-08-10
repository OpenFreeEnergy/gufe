# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import importlib.resources
import urllib.request
from urllib.error import URLError
import functools
import pooch
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from openff.units import unit
import gufe
from gufe.tests.test_protocol import DummyProtocol

try:
    urllib.request.urlopen("https://google.com")
except URLError:
    HAS_INTERNET = False
else:
    HAS_INTERNET = True

PLB_files = pooch.create(
    path=pooch.os_cache('pdbinf'),
    base_url='https://github.com/openforcefield/protein-ligand-benchmark/raw/d3387602bbeb0167abf00dfb81753d8936775dd2/data/',
    version=None,
    registry={
        'p38/01_protein/crd/protein.pdb': '3f0bf718644e7c29f5200cd3def4240ac25ef5fb1948b2e64deb5015d8a45aa4',
        'mcl1/01_protein/crd/protein.pdb': 'f80ff9dd93a5d9dd6e90091e9631a8ce7fe0dc931e16543e22c1f92009660306',
        'cdk2/01_protein/crd/protein.pdb': '15d1e509d7951ca45ea266d51a627d5f452dcf0bb5bd48751ae57eb29e28ab69',
        'shp2/01_protein/crd/protein.pdb': 'd6759cbd135aaddaa658446064df4095d978d3681c014a0528b542d60b2c8770',
        'pde2/01_protein/crd/protein.pdb': '3b7967c1717789215452cdf919520625602d5438a9d2a18620726b8b1b3a8ef0',
        'cmet/01_protein/crd/protein.pdb': '155ec32941a9082dbdbbfde460ff97c88d4fe7e100e9a9577edb5a9e7b6467ae',
        'ptp1b/01_protein/crd/protein.pdb': 'bfa0f9204e96aa463b80946b788c4153cd24701291007eb77638a16fd156634e',
        'thrombin/01_protein/crd/protein.pdb': 'eb4ea18bef9c4c71dcdc922616d6719ee918112be87a0bd6b274c856eff1dd59',
        'cdk8/01_protein/crd/protein.pdb': 'b058774526a19775d8f438b14e9d6da331b6de74e0ef9e96db575f6c0bb067b2',
        'pfkfb3/01_protein/crd/protein.pdb': '4367710db0dbf284cc715ae9a8dd82d06bd77dcc3fb0885678e16632a2732dcc',
        'tyk2/01_protein/crd/protein.pdb': '9090684f4bdae90afbe5f2698a14c778396c024c19ceb6333de4808d9e29fae6',
        'syk/01_protein/crd/protein.pdb': 'f6199d0c1818eb5bb24e164426789cf39cae7aa32c8ca2e98f5f44d299a6f82f',
        'tnks2/01_protein/crd/protein.pdb': 'fc7681a05dbf07590aa8de133f981b6d8ae9cebcc23d54addc2c4fe80be80299',
        'eg5/01_protein/crd/protein.pdb': 'f2964a785c922502dc86fb4e2e5295d32d41d5b68b8c3246e989de5234c3fd0f',
        'hif2a/01_protein/crd/protein.pdb': '5bbf520e7c102a65cc7ba0253fd66f43562f77284c82b3b9613e997b7ac76c93',

    },
)


@pytest.fixture(params=['p38', 'mcl1', 'cdk2', 'shp2', 'pde2', 'cmet', 'ptp1b',
                        'thrombin', 'cdk8', 'pfkfb3', 'tyk2', 'syk', 'tnks2',
                        'eg5', 'hif2a', '181l'])
def PDB_files(request):
    if request.param == '181l':
        with importlib.resources.path('gufe.tests.data', '181l.pdb') as file:
            return str(file)
    else:
        if not HAS_INTERNET:  # pragma: no-cover
            pytest.skip("Skipping because internet seems faulty")
        return PLB_files.fetch('{}/01_protein/crd/protein.pdb'.format(request.param))


@pytest.fixture
def ethane_sdf():
    with importlib.resources.path("gufe.tests.data", "ethane.sdf") as f:
        yield str(f)


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
def offxml_settings_path():
    with importlib.resources.path('gufe.tests.data', 'offxml_settings.json') as f:
        yield str(f)


@pytest.fixture
def all_settings_path():
    with importlib.resources.path('gufe.tests.data', 'all_settings.json') as f:
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


@pytest.fixture(scope='session')
def benzene_modifications():
    with importlib.resources.path('gufe.tests.data',
                                  'benzene_modifications.sdf') as f:
        supp = Chem.SDMolSupplier(str(f), removeHs=False)

        mols = list(supp)

    return {m.GetProp('_Name'): m for m in mols}


@pytest.fixture(scope='session')
def benzene_transforms(benzene_modifications):
    return {k: gufe.SmallMoleculeComponent(v)
            for k, v in benzene_modifications.items()}


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
        protocol=DummyProtocol(settings=None),
        mapping=None,
    )


@pytest.fixture
def complex_equilibrium(solvated_complex):
    return gufe.NonTransformation(
        solvated_complex,
        protocol=DummyProtocol(settings=None)
    )


@pytest.fixture
def benzene_variants_star_map(
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

    variants = [toluene, phenol, benzonitrile, anisole, benzaldehyde,
                styrene]

    # define the solvent chemical systems and transformations between
    # benzene and the others
    solvated_ligands = {}
    solvated_ligand_transformations = {}

    solvated_ligands["benzene"] = gufe.ChemicalSystem(
        {"solvent": solv_comp, "ligand": benzene}, name="benzene-solvent"
    )

    for ligand in variants:
        solvated_ligands[ligand.name] = gufe.ChemicalSystem(
            {"solvent": solv_comp, "ligand": ligand},
            name=f"{ligand.name}-solvnet"
        )
        solvated_ligand_transformations[
            ("benzene", ligand.name)
        ] = gufe.Transformation(
            solvated_ligands["benzene"],
            solvated_ligands[ligand.name],
            protocol=DummyProtocol(settings=None),
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
        solvated_complex_transformations[
            ("benzene", ligand.name)
        ] = gufe.Transformation(
            solvated_complexes["benzene"],
            solvated_complexes[ligand.name],
            protocol=DummyProtocol(settings=None),
            mapping=None,
        )

    return gufe.AlchemicalNetwork(
        list(solvated_ligand_transformations.values())
        + list(solvated_complex_transformations.values())
    )
