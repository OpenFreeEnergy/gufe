# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import os

import pytest
from openff.units import unit

import gufe
from gufe.protocols import ProtocolDAG, execute_DAG
from gufe.protocols.protocolunit import Context
from gufe.storage.externalresource.filestorage import FileStorage
from gufe.storage.storagemanager import StorageManager


class WriterUnit(gufe.ProtocolUnit):
    @staticmethod
    def _execute(ctx: Context, **inputs):
        my_id = inputs["identity"]

        unit_shared_name = f"unit_{my_id}_shared.txt"
        ctx.shared.register(unit_shared_name)
        unit_shared = ctx.scratch / unit_shared_name

        unit_scratch_name = f"unit_{my_id}_scratch.txt"
        unit_scratch = ctx.scratch / unit_scratch_name

        with open(unit_shared, "w") as out:
            out.write(f"unit {my_id} existed\n")
        with open(unit_scratch, "w") as out:
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
def test_execute_dag(tmp_path, keep_shared, keep_scratch, writefile_dag: ProtocolDAG, capture_stderr_stdout):
    scratch = tmp_path / "scratch"
    scratch.mkdir(parents=True)

    shared = tmp_path / "shared"
    shared.mkdir(parents=True)

    shared_storage = FileStorage(shared)

    perm = tmp_path / "perm"
    perm.mkdir(parents=True)
    perm_storage = FileStorage(perm)

    stderr = None
    stdout = None
    if capture_stderr_stdout:
        stderr = tmp_path / "stderr"
        stderr.mkdir(parents=True)
        stdout = tmp_path / "stoud"
        stdout.mkdir(parents=True)

    # run dag
    execute_DAG(
        writefile_dag,
        shared_storage=shared_storage,
        perm_storage=perm_storage,
        scratch_basedir=scratch,
        stderr_basedir=stderr,
        stdout_basedir=stdout,
        keep_shared=keep_shared,
        keep_scratch=keep_scratch,
    )

    # check outputs are as expected
    # will have produced 4 files in scratch and shared directory
    dag_label = str(writefile_dag.key)
    for pu in writefile_dag.protocol_units:
        identity = pu.inputs["identity"]
        # shared_file = os.path.join(shared, f"shared_{str(pu.key)}_attempt_0", f"unit_{identity}_shared.txt")
        shared_file = StorageManager.append_to_namespace(f"{dag_label}/{pu.key}", f"unit_{identity}_shared.txt")
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
            assert shared_storage.exists(shared_file)
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
    dependent_units = [WriterUnit(identity=i, setup=setup_unit, name=f"cycle_{i}") for i in range(1, 4)]

    # Attempt to create a ProtocolDAG without including the setup_unit
    # This should raise a ProtocolDAGError
    with pytest.raises(gufe.protocols.ProtocolDAGError, match="units that were not explicitly provided"):
        gufe.ProtocolDAG(
            protocol_units=dependent_units,  # Missing setup_unit!
            transformation_key=None,
        )
