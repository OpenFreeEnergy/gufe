# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from openff.units import unit
import pytest

import gufe
from .test_tokenization import GufeTokenizableTestsMixin


class DummyProtocolResult(gufe.ProtocolResult):
    def get_estimate(self):
        return self.data['estimate']

    def get_uncertainty(self):
        return self.data['uncertainty']


class TestProtocolResult(GufeTokenizableTestsMixin):
    cls = DummyProtocolResult
    repr = None

    @pytest.fixture
    def instance(self):
        return DummyProtocolResult(
            estimate=4.2 * unit.kilojoule_per_mole,
            uncertainty=0.2 * unit.kilojoule_per_mole,
        )
