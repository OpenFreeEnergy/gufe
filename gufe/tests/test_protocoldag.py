# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import os
import pathlib

import pytest
from openff.units import unit

import gufe
from gufe.protocols import execute_DAG, protocoldag


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
@pytest.mark.parametrize("cache_unitresults", [False, True])
@pytest.mark.parametrize("capture_stderr_stdout", [False, True])
def test_execute_dag(tmp_path, keep_shared, keep_scratch, cache_unitresults, writefile_dag, capture_stderr_stdout):
    shared = pathlib.Path(tmp_path / "shared")
    shared.mkdir(parents=True)

    scratch = pathlib.Path(tmp_path / "scratch")
    scratch.mkdir(parents=True)

    cache_basedir = pathlib.Path(tmp_path / "openfe_cache")
    cache_basedir.mkdir(parents=True)

    stderr = None
    stdout = None
    if capture_stderr_stdout:
        stderr = pathlib.Path(tmp_path / "stderr")
        stderr.mkdir(parents=True)
        stdout = pathlib.Path(tmp_path / "stdout")
        stdout.mkdir(parents=True)

    # run dag
    execute_DAG(
        writefile_dag,
        shared_basedir=shared,
        scratch_basedir=scratch,
        cache_basedir=cache_basedir,
        stderr_basedir=stderr,
        stdout_basedir=stdout,
        keep_shared=keep_shared,
        keep_scratch=keep_scratch,
        cache_unitresults=cache_unitresults,
    )
    # check outputs are as expected
    # will have produced 4 files in scratch and shared directory
    for pu in writefile_dag.protocol_units:
        id = pu.inputs["identity"]
        shared_file = os.path.join(shared, f"shared_{str(pu.key)}_attempt_0", f"unit_{id}_shared.txt")
        scratch_file = os.path.join(scratch, f"scratch_{str(pu.key)}_attempt_0", f"unit_{id}_scratch.txt")
        # TODO: add result key.json
        unit_result_file = os.path.join(cache_basedir, f"unitresults_{str(writefile_dag.key)}")

        if capture_stderr_stdout:
            stderr_file = os.path.join(
                stderr,
                f"stderr_{str(pu.key)}_attempt_0",
                f"unit_{id}_stderr",
            )
            stdout_file = os.path.join(stdout, f"stdout_{str(pu.key)}_attempt_0", f"unit_{id}_stdout")
            # TODO: add result key.json
            unit_result_file = os.path.join(cache_basedir, f"unitresults_{str(writefile_dag.key)}")

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
        if cache_unitresults:
            assert os.path.exists(unit_result_file)
        else:
            assert not os.path.exists(unit_result_file)
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


def test_execute_DAG_cached_unitresults(tmp_path):
    """Test that execute_DAG will re-run based on cache_basedir where only a terminal node is missing results."""

    # Create a setup unit that other units depend on
    setup_unit = WriterUnit(identity=0, name="setup")

    # Create units that depend on the setup unit
    dependent_units = [WriterUnit(identity=i, setup=setup_unit, name=f"cycle_{i}") for i in range(1, 4)]

    dep_dag = gufe.ProtocolDAG(
        protocol_units=dependent_units + [setup_unit],
        transformation_key=None,
    )

    # run all unit_results
    shared = pathlib.Path(tmp_path / "shared")
    shared.mkdir(parents=True)

    scratch = pathlib.Path(tmp_path / "scratch")
    scratch.mkdir(parents=True)

    unit_results_dir = pathlib.Path(tmp_path / "unitresults_cache")
    protocol_result = execute_DAG(
        dep_dag,
        shared_basedir=shared,
        scratch_basedir=scratch,
        cache_basedir=unit_results_dir,
        stderr_basedir=None,
        stdout_basedir=None,
        keep_shared=False,
        keep_scratch=False,
        cache_unitresults=True,
    )

    for pu in dep_dag.protocol_units:
        assert os.path.exists(os.path.join(unit_results_dir, f"unitresults_{dep_dag.key}", f"{str(pu.key)}.json"))

    # choose a terminal result so that only one node is rerun
    pu_to_corrupt = dependent_units[0]

    with open(os.path.join(unit_results_dir, f"unitresults_{dep_dag.key}", f"{str(pu_to_corrupt.key)}.json"), "a") as f:
        f.write("string that will break JSON.")

    protocol_result_rerun = execute_DAG(
        dep_dag,
        shared_basedir=shared,
        scratch_basedir=scratch,
        cache_basedir=unit_results_dir,
        stderr_basedir=None,
        stdout_basedir=None,
        keep_shared=False,
        keep_scratch=False,
        cache_unitresults=True,
    )

    assert protocol_result.protocol_units == protocol_result_rerun.protocol_units
    # if the cache isn't used, these would be identical

    rerun_keys = {r.key for r in protocol_result_rerun.protocol_unit_results}
    original_keys = {r.key for r in protocol_result.protocol_unit_results}

    # Only one result should differ (the corrupted one)
    assert len(rerun_keys.symmetric_difference(original_keys)) == 2
    assert len(rerun_keys.intersection(original_keys)) == len(protocol_result.protocol_unit_results) - 1

    assert protocol_result_rerun.graph.edges == protocol_result.graph.edges


def test_get_valid_unit_results(tmp_path):
    """
    Create a graph of dependencies that looks like this:
    A<-B, B<-C, B<-D, B<-E, D<-F, E<-F
    or read top-down:
        A
        B
      C D E
         F
    """

    unit_A = WriterUnit(identity="A", name="unit_A")
    unit_B = WriterUnit(identity="B", name="unit_B", needs=[unit_A])
    unit_C = WriterUnit(identity="C", name="unit_C", needs=[unit_B])
    unit_D = WriterUnit(identity="D", name="unit_D", needs=[unit_B])
    unit_E = WriterUnit(identity="E", name="unit_E", needs=[unit_B])
    unit_F = WriterUnit(identity="F", name="unit_F", needs=[unit_D, unit_E])

    all_protocol_units = {unit_A, unit_B, unit_C, unit_D, unit_E, unit_F}
    # Create units that depend on the setup unit
    dep_dag = gufe.ProtocolDAG(
        protocol_units=all_protocol_units,
        transformation_key=None,
    )
    shared = pathlib.Path(tmp_path / "shared")
    shared.mkdir(parents=True)

    scratch = pathlib.Path(tmp_path / "scratch")
    scratch.mkdir(parents=True)

    unit_results_dir = pathlib.Path(tmp_path / "unitresults_cache")
    protocol_result = execute_DAG(
        dep_dag,
        shared_basedir=shared,
        scratch_basedir=scratch,
        cache_basedir=unit_results_dir,
        stderr_basedir=None,
        stdout_basedir=None,
        keep_shared=False,
        keep_scratch=False,
        cache_unitresults=True,
    )
    all_cached_unit_results = protocol_result.protocol_unit_results

    # cache is empty, so nothing should be skipped
    units_to_skip = protocoldag._get_valid_unit_results(protocoldag=dep_dag, unit_results=[])
    assert len(units_to_skip.keys()) == 0

    # all results are available, so everything should be skipped
    units_to_skip = protocoldag._get_valid_unit_results(protocoldag=dep_dag, unit_results=all_cached_unit_results)
    assert set(units_to_skip.keys()) == {u.key for u in all_protocol_units}

    # drop the top-most unit, so nothing should be skipped
    unit_results_drop_A = [u for u in all_cached_unit_results if u.name != "unit_A"]
    units_to_skip = protocoldag._get_valid_unit_results(protocoldag=dep_dag, unit_results=unit_results_drop_A)
    assert len(units_to_skip.keys()) == 0

    # drop terminal nodes, so everything *but* the terminal nodes can be skipped
    unit_results_drop_C_F = [u for u in all_cached_unit_results if u.name not in ("unit_C", "unit_F")]
    units_to_skip = protocoldag._get_valid_unit_results(protocoldag=dep_dag, unit_results=unit_results_drop_C_F)
    assert set(units_to_skip.keys()) == {u.key for u in (unit_A, unit_B, unit_D, unit_E)}
