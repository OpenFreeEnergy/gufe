# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import pytest

from gufe.transformations import Transformation, NonTransformation
from gufe.protocols.protocoldag import execute_DAG

from .test_protocol import DummyProtocol, DummyProtocolResult


@pytest.fixture
def absolute_transformation(solvated_ligand, solvated_complex):
    return Transformation(
        solvated_ligand,
        solvated_complex,
        protocol=DummyProtocol(settings=None),
        mapping=None,
    )


@pytest.fixture
def complex_equilibrium(solvated_complex):
    return NonTransformation(solvated_complex, protocol=DummyProtocol(settings=None))


class TestTransformation:
    def test_init(self, absolute_transformation, solvated_ligand, solvated_complex):
        tnf = absolute_transformation

        assert tnf.stateA is solvated_ligand
        assert tnf.stateB is solvated_complex

    def test_protocol(self, absolute_transformation):
        tnf = absolute_transformation

        assert isinstance(tnf.protocol, DummyProtocol)

        protocoldag = tnf.create()
        protocoldagresult = execute_DAG(protocoldag)

        protocolresult = tnf.gather([protocoldagresult])

        assert isinstance(protocolresult, DummyProtocolResult)

        len(protocolresult.data) == 1

    def test_equality(self, absolute_transformation, solvated_ligand, solvated_complex):

        opposite = Transformation(
            solvated_complex, solvated_ligand, protocol=DummyProtocol(settings=None)
        )
        assert absolute_transformation != opposite

        different_protocol_settings = Transformation(
            solvated_ligand,
            solvated_complex,
            protocol=DummyProtocol(settings={"lol": True}),
        )
        assert absolute_transformation != different_protocol_settings

        identical = Transformation(
            solvated_ligand,
            solvated_complex,
            protocol=DummyProtocol(settings=None),
            mapping=None,
        )
        assert absolute_transformation == identical

    def test_dict_roundtrip(self):
        # TODO: need registration of `Protocol`s for this to work
        ...


class TestNonTransformation:
    def test_init(self, complex_equilibrium, solvated_complex):

        ntnf = complex_equilibrium

        assert ntnf.system is solvated_complex

    def test_protocol(self, complex_equilibrium):
        ntnf = complex_equilibrium

        assert isinstance(ntnf.protocol, DummyProtocol)

        protocoldag = ntnf.create()
        protocoldagresult = execute_DAG(protocoldag)

        protocolresult = ntnf.gather([protocoldagresult])

        assert isinstance(protocolresult, DummyProtocolResult)

        len(protocolresult.data) == 1

    def test_equality(self, complex_equilibrium, solvated_ligand, solvated_complex):

        different_protocol_settings = NonTransformation(
            solvated_complex, protocol=DummyProtocol(settings={"lol": True})
        )
        assert complex_equilibrium != different_protocol_settings

        identical = NonTransformation(
            solvated_complex, protocol=DummyProtocol(settings=None)
        )
        assert complex_equilibrium == identical

        different_system = NonTransformation(
            solvated_ligand, protocol=DummyProtocol(settings=None)
        )
        assert complex_equilibrium != different_system

    def test_dict_roundtrip(self):
        # TODO: need registration of `Protocol`s for this to work
        ...
