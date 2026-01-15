# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import os
import pathlib

import pytest
from openff.units import unit

import gufe
from gufe.protocols import execute_DAG


class WriterUnit(gufe.ProtocolUnit):
    @staticmethod
    def _execute(ctx, **inputs):
        my_id = inputs["identity"]

        with open(os.path.join(ctx.shared, f"unit_{my_id}_shared.txt"), "w") as out:
            out.write(f"unit {my_id} existed!\n")
        with open(os.path.join(ctx.scratch, f"unit_{my_id}_scratch.txt"), "w") as out:
            out.write(f"unit {my_id} was here\n")

        if ctx.stderr:
            with open(os.path.join(ctx.stderr, f"unit_{my_id}_stderr"), "w") as out:
                out.write(f"unit {my_id} wrote to stderr")
        if ctx.stdout:
            with open(os.path.join(ctx.stdout, f"unit_{my_id}_stdout"), "w") as out:
                out.write(f"unit {my_id} wrote to stdout")

        return {
            "log": "finished",
        }


class WriterSettings(gufe.settings.Settings):
    n_repeats: int


class WriterProtocolResult(gufe.ProtocolResult):
    def get_estimate(self): ...

    def get_uncertainty(self): ...


class WriterProtocol(gufe.Protocol):
    result_cls = WriterProtocolResult
    _settings_cls = WriterSettings

    @classmethod
    def _default_settings(cls):
        return WriterSettings(
            thermo_settings=gufe.settings.ThermoSettings(temperature=298 * unit.kelvin),
            forcefield_settings=gufe.settings.OpenMMSystemGeneratorFFSettings(),
            n_repeats=4,
        )

    @classmethod
    def _defaults(cls):
        return {}

    def _create(self, stateA, stateB, mapping, extends=None) -> list[gufe.ProtocolUnit]:
        return [WriterUnit(identity=i) for i in range(self.settings.n_repeats)]  # type: ignore

    def _gather(self, results):
        return {}


@pytest.fixture()
def writefile_dag():
    s1 = gufe.ChemicalSystem(components={})
    s2 = gufe.ChemicalSystem(components={})

    p = WriterProtocol(settings=WriterProtocol.default_settings())

    return p.create(stateA=s1, stateB=s2, mapping=[])


@pytest.mark.parametrize("keep_shared", [False, True])
@pytest.mark.parametrize("keep_scratch", [False, True])
@pytest.mark.parametrize("capture_stderr_stdout", [False, True])
def test_execute_dag(tmpdir, keep_shared, keep_scratch, writefile_dag, capture_stderr_stdout):
    with tmpdir.as_cwd():
        shared = pathlib.Path("shared")
        shared.mkdir(parents=True)

        scratch = pathlib.Path("scratch")
        scratch.mkdir(parents=True)

        stderr = None
        stdout = None
        if capture_stderr_stdout:
            stderr = pathlib.Path("stderr")
            stderr.mkdir(parents=True)
            stdout = pathlib.Path("stdout")
            stdout.mkdir(parents=True)

        # run dag
        execute_DAG(
            writefile_dag,
            shared_basedir=shared,
            scratch_basedir=scratch,
            stderr_basedir=stderr,
            stdout_basedir=stdout,
            keep_shared=keep_shared,
            keep_scratch=keep_scratch,
        )

        # check outputs are as expected
        # will have produced 4 files in scratch and shared directory
        for pu in writefile_dag.protocol_units:
            identity = pu.inputs["identity"]
            shared_file = os.path.join(shared, f"shared_{str(pu.key)}_attempt_0", f"unit_{identity}_shared.txt")
            scratch_file = os.path.join(
                scratch,
                f"scratch_{str(pu.key)}_attempt_0",
                f"unit_{identity}_scratch.txt",
            )

            if capture_stderr_stdout:
                stderr_file = os.path.join(
                    stderr,
                    f"stderr_{str(pu.key)}_attempt_0",
                    f"unit_{identity}_stderr",
                )
                stdout_file = os.path.join(
                    stdout,
                    f"stdout_{str(pu.key)}_attempt_0",
                    f"unit_{identity}_stdout",
                )

                # stderr and stdout are always removed since their
                # contents are included in the unit results
                assert not os.path.exists(stderr_file)
                assert not os.path.exists(stdout_file)

            if keep_shared:
                assert os.path.exists(shared_file)
            else:
                assert not os.path.exists(shared_file)
            if keep_scratch:
                assert os.path.exists(scratch_file)
            else:
                assert not os.path.exists(scratch_file)

        # check that our shared and scratch basedirs are left behind
        assert shared.exists()
        assert scratch.exists()


def test_protocoldag_missing_dependency_unit():
    """Test that ProtocolDAG raises an error when units with dependencies
    are provided but their dependencies are not explicitly included.

    This test addresses issue #583: Protocol._create should return all units,
    not rely on implicit dependency discovery.
    """
    # Create a setup unit that other units depend on
    setup_unit = WriterUnit(identity=0, name="setup")

    # Create units that depend on the setup unit
    dependent_units = [
        WriterUnit(identity=i, setup=setup_unit, name=f"cycle_{i}")
        for i in range(1, 4)
    ]

    # Attempt to create a ProtocolDAG without including the setup_unit
    # This should raise a ProtocolDAGError
    with pytest.raises(gufe.protocols.ProtocolDAGError, match="units that were not explicitly provided"):
        gufe.ProtocolDAG(
            protocol_units=dependent_units,  # Missing setup_unit!
            transformation_key=None,
        )
