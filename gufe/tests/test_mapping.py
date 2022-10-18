# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import gufe
from gufe import AtomMapping
import pytest


from .test_tokenization import GufeTokenizableTestsMixin


class ExampleMapping(AtomMapping):
    def __init__(self, molA: gufe.SmallMoleculeComponent,
                 molB: gufe.SmallMoleculeComponent, mapping):
        super().__init__(molA, molB)
        self._mapping = mapping

    def _defaults(self):
        return {}

    def _to_dict(self):
        return {
            'molA': self._componentA,
            'molB': self._componentB,
            'mapping': self._mapping,
        }

    @classmethod
    def _from_dict(cls, d):
        return cls(**d)

    def componentA_to_componentB(self):
        return self._mapping

    def componentB_to_componentA(self):
        return {v: k for k, v in self._mapping}

    def componentA_unique(self):
        return (i for i in range(self._molA.to_rdkit().GetNumAtoms())
                if i not in self._mapping)

    def componentB_unique(self):
        return (i for i in range(self._molB.to_rdkit().GetNumAtoms())
                if i not in self._mapping.values())


class TestMappingAbstractClass(GufeTokenizableTestsMixin):
    cls = ExampleMapping
    key = 'ExampleMapping-43716e2d20e7c33b99ea470707c5832f'

    @pytest.fixture
    def instance(self, benzene, toluene):
        indices = {1: 5, 2: 6}

        return ExampleMapping(benzene, toluene, indices)

    def test_contains(self, instance, benzene, toluene, phenol):
        assert benzene in instance
        assert toluene in instance
        assert phenol not in instance

    def test_component_access(self, instance, benzene, toluene):
        assert instance.componentA == benzene
        assert instance.componentB == toluene
