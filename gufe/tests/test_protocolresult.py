# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from openff.units import unit
import pytest

import gufe
from .test_tokenization import GufeTokenizableTestsMixin


class DummyProtocolResult(gufe.ProtocolResult):
    def get_estimate(self):
        return self.data["estimate"]

    def get_uncertainty(self):
        return self.data["uncertainty"]


class TestProtocolResult(GufeTokenizableTestsMixin):
    cls = DummyProtocolResult
    key = "DummyProtocolResult-38c5649090dae613b2570c28155cc5fa"
    repr = f"<{key}>"

    @pytest.fixture
    def instance(self):
        return DummyProtocolResult(
            estimate=4.2 * unit.kilojoule_per_mole,
            uncertainty=0.2 * unit.kilojoule_per_mole,
        )

    def test_protocolresult_get_estimate(self, instance):
        assert instance.get_estimate() == 4.2 * unit.kilojoule_per_mole

    def test_protocolresult_get_uncertainty(self, instance):
        assert instance.get_uncertainty() == 0.2 * unit.kilojoule_per_mole

    def test_protocolresult_default_n_protocol_dag_results(self, instance):
        assert instance.n_protocol_dag_results == 0

    @pytest.mark.parametrize(
        "arg, expected", [(0, 0), (1, 1), (-1, ValueError)]
    )
    def test_protocolresult_get_n_protocol_dag_results_args(self, arg, expected):
        try:
            protocol_result = DummyProtocolResult(
                n_protocol_dag_results=arg,
                estimate=4.2 * unit.kilojoule_per_mole,
                uncertainty=0.2 * unit.kilojoule_per_mole,
            )
            assert protocol_result.n_protocol_dag_results == expected
        except ValueError:
            if expected is not ValueError:
                raise AssertionError()
