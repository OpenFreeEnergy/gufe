# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import pytest
from gufe import AtomMappingScorer, AtomMapper


def test_atom_mapping_scorer():
    with pytest.raises(TypeError, match="Can't instantiate abstract class AtomMappingScorer without an implementation for abstract methods '_defaults', '_from_dict', '_to_dict', 'get_score'"):
        scorer = AtomMappingScorer()

def test_atom_mapper():
    with pytest.raises(TypeError, match="Can't instantiate abstract class AtomMapper without an implementation for abstract methods '_defaults', '_from_dict', '_to_dict', 'suggest_mappings'"):
        mapper = AtomMapper()