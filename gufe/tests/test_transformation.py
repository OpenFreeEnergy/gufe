# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import io
import pathlib

import pytest

import gufe
from gufe.protocols.errors import ProtocolValidationError
from gufe.protocols.protocoldag import execute_DAG
from gufe.transformations import NonTransformation, Transformation

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
    return NonTransformation(
        solvated_complex,
        protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
    )


class TestTransformation(GufeTokenizableTestsMixin):
    cls = Transformation
    repr = "Transformation(stateA=ChemicalSystem(name=, components={'ligand': SmallMoleculeComponent(name=toluene), 'solvent': SolventComponent(name=O, K+, Cl-)}), stateB=ChemicalSystem(name=, components={'protein': ProteinComponent(name=), 'solvent': SolventComponent(name=O, K+, Cl-), 'ligand': SmallMoleculeComponent(name=toluene)}), protocol=<DummyProtocol-13cfc31d879fd8877e6a57558c2e9f38>, name=None)"

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
            shared = pathlib.Path("shared")
            shared.mkdir(parents=True)

            scratch = pathlib.Path("scratch")
            scratch.mkdir(parents=True)

            stderr = pathlib.Path("stderr")
            stderr.mkdir(parents=True)

            stdout = pathlib.Path("stdout")
            stdout.mkdir(parents=True)

            protocoldagresult = execute_DAG(
                protocoldag,
                shared_basedir=shared,
                scratch_basedir=scratch,
                stderr_basedir=stderr,
                stdout_basedir=stdout,
            )

        protocolresult = tnf.gather([protocoldagresult])

        assert isinstance(protocolresult, DummyProtocolResult)

        assert len(protocolresult.data) == 2
        assert "logs" in protocolresult.data
        assert "key_results" in protocolresult.data

    def test_validation_good_inputs(self, solvated_ligand, solvated_complex):
        Transformation(
            solvated_ligand,
            solvated_complex,
            protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
            mapping=None,
            validate=True,
        )

    def test_validation_bad_state(self, solvated_ligand, solvated_complex, tmpdir):
        # show that instantiation without validation is possible (even though it's wrong)
        Transformation(
            solvated_ligand,
            None,
            protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
            mapping=None,
        )

        with pytest.raises(
            ProtocolValidationError, match="`stateA` and `stateB` must be instances of a `ChemicalSystem`"
        ):
            tf = Transformation(
                solvated_ligand,
                None,
                protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
                mapping=None,
                validate=True,
            )

    def test_protocol_extend(self, absolute_transformation, tmpdir):
        tnf = absolute_transformation

        assert isinstance(tnf.protocol, DummyProtocol)

        with tmpdir.as_cwd():
            shared = pathlib.Path("shared")
            shared.mkdir(parents=True)

            scratch = pathlib.Path("scratch")
            scratch.mkdir(parents=True)

            stderr = pathlib.Path("stderr")
            stderr.mkdir(parents=True)

            stdout = pathlib.Path("stdout")
            stdout.mkdir(parents=True)

            protocoldag = tnf.create()
            protocoldagresult = execute_DAG(
                protocoldag,
                shared_basedir=shared,
                scratch_basedir=scratch,
                stderr_basedir=stderr,
                stdout_basedir=stdout,
            )

            protocoldag2 = tnf.create(extends=protocoldagresult)
            protocoldagresult2 = execute_DAG(
                protocoldag2,
                shared_basedir=shared,
                scratch_basedir=scratch,
                stderr_basedir=stderr,
                stdout_basedir=stdout,
            )

        protocolresult = tnf.gather([protocoldagresult, protocoldagresult2])

        assert isinstance(protocolresult, DummyProtocolResult)

        assert len(protocolresult.data) == 2

    def test_equality(self, absolute_transformation, solvated_ligand, solvated_complex):
        opposite = Transformation(
            solvated_complex,
            solvated_ligand,
            protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
        )
        assert absolute_transformation != opposite

        s = DummyProtocol.default_settings()
        s.n_repeats = 99
        different_protocol_settings = Transformation(
            solvated_ligand,
            solvated_complex,
            protocol=DummyProtocol(settings=s),
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
        with pytest.warns(DeprecationWarning, match="use of this method is deprecated; instead use `to_json`"):
            absolute_transformation.dump(string)
        string.seek(0)
        with pytest.warns(DeprecationWarning, match="use of this method is deprecated; instead use `from_json`"):
            recreated = Transformation.load(string)
        assert absolute_transformation == recreated

    def test_deprecation_warning_on_dict_mapping(self, solvated_ligand, solvated_complex):
        lig = solvated_complex.components["ligand"]
        # this mapping makes no sense, but it'll trigger the dep warning we want
        mapping = gufe.LigandAtomMapping(lig, lig, componentA_to_componentB={})

        with pytest.warns(DeprecationWarning, match="mapping input as a dict is deprecated"):
            Transformation(
                solvated_complex,
                solvated_ligand,
                protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
                mapping={"ligand": mapping},
            )


class TestNonTransformation(GufeTokenizableTestsMixin):
    cls = NonTransformation
    repr = "NonTransformation(system=ChemicalSystem(name=, components={'protein': ProteinComponent(name=), 'solvent': SolventComponent(name=O, K+, Cl-), 'ligand': SmallMoleculeComponent(name=toluene)}), protocol=<DummyProtocol-13cfc31d879fd8877e6a57558c2e9f38>, name=None)"

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
            shared = pathlib.Path("shared")
            shared.mkdir(parents=True)

            scratch = pathlib.Path("scratch")
            scratch.mkdir(parents=True)

            stderr = pathlib.Path("stderr")
            stderr.mkdir(parents=True)

            stdout = pathlib.Path("stdout")
            stdout.mkdir(parents=True)

            protocoldagresult = execute_DAG(
                protocoldag,
                shared_basedir=shared,
                scratch_basedir=scratch,
                stderr_basedir=stderr,
                stdout_basedir=stdout,
            )

        protocolresult = ntnf.gather([protocoldagresult])

        assert isinstance(protocolresult, DummyProtocolResult)

        assert len(protocolresult.data) == 2
        assert "logs" in protocolresult.data
        assert "key_results" in protocolresult.data

    def test_validation_good_inputs(self, solvated_complex):
        NonTransformation(
            solvated_complex,
            protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
        )

    def test_validation_bad_state(self, solvated_ligand, solvated_complex, tmpdir):
        # show that instantiation without validation is possible (even though it's wrong)
        NonTransformation(
            None,
            protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
        )

        with pytest.raises(
            ProtocolValidationError, match="`stateA` and `stateB` must be instances of a `ChemicalSystem`"
        ):
            ntf = NonTransformation(
                None,
                protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
                validate=True,
            )

    def test_protocol_extend(self, complex_equilibrium, tmpdir):
        ntnf = complex_equilibrium

        assert isinstance(ntnf.protocol, DummyProtocol)

        with tmpdir.as_cwd():
            shared = pathlib.Path("shared")
            shared.mkdir(parents=True)

            scratch = pathlib.Path("scratch")
            scratch.mkdir(parents=True)

            stderr = pathlib.Path("stderr")
            stderr.mkdir(parents=True)

            stdout = pathlib.Path("stdout")
            stdout.mkdir(parents=True)

            protocoldag = ntnf.create()
            protocoldagresult = execute_DAG(
                protocoldag,
                shared_basedir=shared,
                scratch_basedir=scratch,
                stderr_basedir=stderr,
                stdout_basedir=stdout,
            )

            protocoldag2 = ntnf.create(extends=protocoldagresult)
            protocoldagresult2 = execute_DAG(
                protocoldag2,
                shared_basedir=shared,
                scratch_basedir=scratch,
            )

        protocolresult = ntnf.gather([protocoldagresult, protocoldagresult2])

        assert isinstance(protocolresult, DummyProtocolResult)

        assert len(protocolresult.data) == 2

    def test_equality(self, complex_equilibrium, solvated_ligand, solvated_complex):
        s = DummyProtocol.default_settings()
        s.n_repeats = 4031
        different_protocol_settings = NonTransformation(solvated_complex, protocol=DummyProtocol(settings=s))
        assert complex_equilibrium != different_protocol_settings

        identical = NonTransformation(
            solvated_complex,
            protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
        )
        assert complex_equilibrium == identical

        different_system = NonTransformation(
            solvated_ligand,
            protocol=DummyProtocol(settings=DummyProtocol.default_settings()),
        )
        assert complex_equilibrium != different_system
