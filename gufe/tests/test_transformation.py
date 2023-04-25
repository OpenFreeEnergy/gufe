# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import pytest
import io
import pathlib

from gufe.transformations import Transformation, NonTransformation
from gufe.protocols.protocoldag import execute_DAG

from .test_protocol import DummyProtocol, DummyProtocolResult
from .test_tokenization import GufeTokenizableTestsMixin


@pytest.fixture
def absolute_transformation(solvated_ligand, solvated_complex):
    return Transformation(
        solvated_ligand,
        solvated_complex,
        protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
        mapping=None,
    )


@pytest.fixture
def complex_equilibrium(solvated_complex):
    return NonTransformation(solvated_complex,
                             protocol=DummyProtocol(settings=DummyProtocol.default_settings()))


class TestTransformation(GufeTokenizableTestsMixin):

    cls = Transformation
    key = "Transformation-e710da52fcf4b7d78f890e6e00e8f6c5"
    repr = "Transformation(stateA=ChemicalSystem(name=, components={'ligand': SmallMoleculeComponent(name=toluene), 'solvent': SolventComponent(name=O, K+, Cl-)}), stateB=ChemicalSystem(name=, components={'protein': ProteinComponent(name=), 'solvent': SolventComponent(name=O, K+, Cl-), 'ligand': SmallMoleculeComponent(name=toluene)}), protocol=<DummyProtocol-f4f0aab0421240c9bc65b6491791f9d8>)"

    @pytest.fixture
    def instance(self, absolute_transformation):
        return absolute_transformation

    def test_init(self, absolute_transformation, solvated_ligand, solvated_complex):
        tnf = absolute_transformation

        assert tnf.stateA is solvated_ligand
        assert tnf.stateB is solvated_complex

    def test_protocol(self, absolute_transformation, tmpdir):
        tnf = absolute_transformation

        assert isinstance(tnf.protocol, DummyProtocol)

        protocoldag = tnf.create()

        with tmpdir.as_cwd():
            shared = pathlib.Path('shared')
            shared.mkdir(parents=True)

            scratch = pathlib.Path('scratch')
            scratch.mkdir(parents=True)
    
            protocoldagresult = execute_DAG(protocoldag, shared_basedir=shared, scratch_basedir=scratch)

        protocolresult = tnf.gather([protocoldagresult])

        assert isinstance(protocolresult, DummyProtocolResult)

        assert len(protocolresult.data) == 2
        assert 'logs' in protocolresult.data
        assert 'key_results' in protocolresult.data

    def test_protocol_extend(self, absolute_transformation, tmpdir):
        tnf = absolute_transformation

        assert isinstance(tnf.protocol, DummyProtocol)

        with tmpdir.as_cwd():
            shared = pathlib.Path('shared')
            shared.mkdir(parents=True)

            scratch = pathlib.Path('scratch')
            scratch.mkdir(parents=True)
    
            protocoldag = tnf.create()
            protocoldagresult = execute_DAG(protocoldag, shared_basedir=shared, scratch_basedir=scratch)

            protocoldag2 = tnf.create(extends=protocoldagresult)
            protocoldagresult2 = execute_DAG(protocoldag2, shared_basedir=shared, scratch_basedir=scratch)

        protocolresult = tnf.gather([protocoldagresult, protocoldagresult2])

        assert isinstance(protocolresult, DummyProtocolResult)

        assert len(protocolresult.data) == 2

    def test_equality(self, absolute_transformation, solvated_ligand, solvated_complex):

        opposite = Transformation(
            solvated_complex, solvated_ligand,
            protocol=DummyProtocol(settings=DummyProtocol.default_settings())
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
            protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
            mapping=None,
        )
        assert absolute_transformation == identical

    def test_dump_load_roundtrip(self, absolute_transformation):
        string = io.StringIO()
        absolute_transformation.dump(string)
        string.seek(0)
        recreated = Transformation.load(string)
        assert absolute_transformation == recreated


class TestNonTransformation(GufeTokenizableTestsMixin):

    cls = NonTransformation
    key = "NonTransformation-55e8c1e7274d71e246d36ca4d718fed3"
    repr = "NonTransformation(stateA=ChemicalSystem(name=, components={'protein': ProteinComponent(name=), 'solvent': SolventComponent(name=O, K+, Cl-), 'ligand': SmallMoleculeComponent(name=toluene)}), stateB=ChemicalSystem(name=, components={'protein': ProteinComponent(name=), 'solvent': SolventComponent(name=O, K+, Cl-), 'ligand': SmallMoleculeComponent(name=toluene)}), protocol=<DummyProtocol-f4f0aab0421240c9bc65b6491791f9d8>)"

    @pytest.fixture
    def instance(self, complex_equilibrium):
        return complex_equilibrium

    def test_init(self, complex_equilibrium, solvated_complex):

        ntnf = complex_equilibrium

        assert ntnf.system is solvated_complex

    def test_protocol(self, complex_equilibrium, tmpdir):
        ntnf = complex_equilibrium

        assert isinstance(ntnf.protocol, DummyProtocol)

        protocoldag = ntnf.create()

        with tmpdir.as_cwd():
            shared = pathlib.Path('shared')
            shared.mkdir(parents=True)

            scratch = pathlib.Path('scratch')
            scratch.mkdir(parents=True)

            protocoldagresult = execute_DAG(protocoldag, shared_basedir=shared, scratch_basedir=scratch)

        protocolresult = ntnf.gather([protocoldagresult])

        assert isinstance(protocolresult, DummyProtocolResult)

        assert len(protocolresult.data) == 2
        assert 'logs' in protocolresult.data
        assert 'key_results' in protocolresult.data

    def test_protocol_extend(self, complex_equilibrium, tmpdir):
        ntnf = complex_equilibrium

        assert isinstance(ntnf.protocol, DummyProtocol)

        with tmpdir.as_cwd():
            shared = pathlib.Path('shared')
            shared.mkdir(parents=True)

            scratch = pathlib.Path('scratch')
            scratch.mkdir(parents=True)

            protocoldag = ntnf.create()
            protocoldagresult = execute_DAG(protocoldag, shared_basedir=shared, scratch_basedir=scratch)

            protocoldag2 = ntnf.create(extends=protocoldagresult)
            protocoldagresult2 = execute_DAG(protocoldag2, shared_basedir=shared, scratch_basedir=scratch)

        protocolresult = ntnf.gather([protocoldagresult, protocoldagresult2])

        assert isinstance(protocolresult, DummyProtocolResult)

        assert len(protocolresult.data) == 2

    def test_equality(self, complex_equilibrium, solvated_ligand, solvated_complex):

        different_protocol_settings = NonTransformation(
            solvated_complex, protocol=DummyProtocol(settings={"lol": True})
        )
        assert complex_equilibrium != different_protocol_settings

        identical = NonTransformation(
            solvated_complex, protocol=DummyProtocol(settings=DummyProtocol.default_settings())
        )
        assert complex_equilibrium == identical

        different_system = NonTransformation(
            solvated_ligand, protocol=DummyProtocol(settings=DummyProtocol.default_settings())
        )
        assert complex_equilibrium != different_system
