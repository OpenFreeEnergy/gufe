# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import importlib.resources
import pytest
from rdkit import Chem

import gufe

@pytest.fixture
def serialization_template():
    def inner(filename):
        loc = "gufe.tests.data"
        tmpl = importlib.resources.read_text(loc, filename)
        return tmpl.format(OFE_VERSION=gufe.__version__)

    return inner

@pytest.fixture(scope='session')
def ethane():
    return gufe.LigandComponent(Chem.MolFromSmiles('CC'))


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
def PDBx_181L_path():
    with importlib.resources.path('gufe.tests.data', '181l.cif') as f:
        yield str(f)
