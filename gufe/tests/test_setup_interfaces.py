# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import pytest
from gufe import AtomMappingScorer, AtomMapper


def test_atom_mapping_scorer():
    with pytest.raises(TypeError, match=""Can't instantiate abstract class AtomMappingScorer"):
        scorer = AtomMappingScorer()


def test_atom_mapper():
    with pytest.raises(TypeError, match=""Can't instantiate abstract class AtomMapper"):
        mapper = AtomMapper()
