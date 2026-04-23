import pathlib
from pathlib import Path
from unittest.mock import patch

from openfe.protocols import openmm_rfe
from openfe.protocols.openmm_rfe import _rfe_utils
from openff.units import unit
from openmmtools import multistate

from gufe.protocols import execute_DAG


class MockUUID:
    """A mock UUID object that mimics uuid.UUID4 behavior"""

    def __init__(self, value: int):
        self.value = value

    @property
    def hex(self) -> str:
        """Return hex representation like a real UUID"""
        return f"hex-{self.value}"

    def __int__(self) -> int:
        """Allow int() conversion for backward compatibility"""
        return self.value

    def __str__(self) -> str:
        """String representation"""
        return f"MockUUID-{self.value}"


class MonotonicUUID:
    """A mock UUID generator that returns monotonically increasing integers"""

    def __init__(self, start=0):
        self.counter = start

    def __call__(self) -> MockUUID:
        """Return a MockUUID object with current counter value and increment"""
        current = MockUUID(self.counter)
        self.counter += 1
        return current

    def reset(self, start=0):
        """Reset the counter to a specific value"""
        self.counter = start


def create_protocol_settings(*, production_length) -> openmm_rfe.RelativeHybridTopologyProtocolSettings:
    settings = openmm_rfe.RelativeHybridTopologyProtocol.default_settings()
    settings.simulation_settings.minimization_steps = 5000  # default is 5000
    settings.simulation_settings.equilibration_length = 10 * unit.picoseconds
    settings.simulation_settings.production_length = production_length
    settings.output_settings.checkpoint_interval = 1 * unit.picoseconds
    settings.protocol_repeats = 1
    return settings


def restore_sampler_from_checkpoint(
    shared_basepath: Path,
) -> _rfe_utils.multistate.HybridRepexSampler:
    nc = shared_basepath / "simulation.nc"
    chk = shared_basepath / "checkpoint.chk"
    reporter = multistate.MultiStateReporter(str(nc), checkpoint_storage=str(chk), open_mode="r")
    sampler = _rfe_utils.multistate.HybridRepexSampler.from_storage(reporter)
    return sampler


def create_protocol_dag(*, benzene_system, toluene_system, benzene_to_toluene_mapping, production_length):
    settings = create_protocol_settings(production_length=production_length)

    protocol = openmm_rfe.RelativeHybridTopologyProtocol(
        settings=settings,
    )

    dag = protocol.create(
        stateA=benzene_system,
        stateB=toluene_system,
        mapping=benzene_to_toluene_mapping,
    )
    return dag


@patch("openfe.protocols.openmm_rfe.equil_rfe_methods.uuid.uuid4")
def test_protocol_dag_resume(mock_uuid, tmpdir, benzene_system, toluene_system, benzene_to_toluene_mapping):
    with tmpdir.as_cwd():
        # Start from 1000
        uuid_generator = MonotonicUUID(start=1000)
        mock_uuid.side_effect = uuid_generator

        shared = pathlib.Path("shared")
        shared.mkdir(parents=True, exist_ok=True)

        scratch = pathlib.Path("scratch")
        scratch.mkdir(parents=True, exist_ok=True)

        protocol_dag = create_protocol_dag(
            benzene_system=benzene_system,
            toluene_system=toluene_system,
            benzene_to_toluene_mapping=benzene_to_toluene_mapping,
            production_length=1 * unit.picoseconds,
        )

        execute_DAG(
            protocol_dag,
            shared_basedir=shared,
            scratch_basedir=scratch,
            keep_shared=True,
            keep_scratch=True,
        )

        # 1 directory in shared
        assert len(list(shared.glob("*"))) == 1

        # Reset the counter, run again
        uuid_generator = MonotonicUUID(start=1000)
        mock_uuid.side_effect = uuid_generator

        execute_DAG(
            protocol_dag,
            shared_basedir=shared,
            scratch_basedir=scratch,
            keep_shared=True,
            keep_scratch=True,
        )

        # 1 directory in shared
        assert len(list(shared.glob("*"))) == 1
