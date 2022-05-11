# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import pytest

import gufe


@pytest.fixture
def absolute_transformation(solvated_ligand, solvated_complex):
    return gufe.Transformation(
            solvated_ligand, 
            solvated_complex, 
            mapping=None,
            protocol=None)


class TestTransformation:
    ...


class TestNonTransformation:
    ...
