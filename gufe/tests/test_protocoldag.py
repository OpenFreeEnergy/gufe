# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import os
import pathlib
import pytest
from openff.units import unit

import gufe
from gufe.protocols import execute_DAG
from gufe.protocols.protocoldag import ReproduceOldBehaviorStorageManager


class WriterUnit(gufe.ProtocolUnit):
    @staticmethod
    def _execute(ctx, **inputs):
        my_id = inputs['identity']

        with open(ctx.shared / f'unit_{my_id}_shared.txt', 'w') as out:
            out.write(f'unit {my_id} existed!\n')
        with open(ctx.scratch / f'unit_{my_id}_scratch.txt', 'w') as out:
            out.write(f'unit {my_id} was here\n')

        return {
            'log': 'finished',
        }


class WriterSettings(gufe.settings.Settings):
    n_repeats: int


class WriterProtocolResult(gufe.ProtocolResult):

    def get_estimate(self):
        ...

    def get_uncertainty(self):
        ...


class WriterProtocol(gufe.Protocol):
    result_cls = WriterProtocolResult

    @classmethod
    def _default_settings(cls):
        return WriterSettings(
            thermo_settings=gufe.settings.ThermoSettings(temperature=298 * unit.kelvin),
            forcefield_settings=gufe.settings.OpenMMSystemGeneratorFFSettings(),
            n_repeats=4
        )
        

    @classmethod
    def _defaults(cls):
        return {}

    def _create(self, stateA,  stateB, mapping=None, extends=None) -> list[gufe.ProtocolUnit]:
        return [
            WriterUnit(identity=i) for i in range(self.settings.n_repeats) # type: ignore
        ]

    def _gather(self, results):
        return {}


@pytest.fixture()
def writefile_dag():
    s1 = gufe.ChemicalSystem(components={})
    s2 = gufe.ChemicalSystem(components={})

    p = WriterProtocol(settings=WriterProtocol.default_settings())

    return p.create(stateA=s1, stateB=s2, mapping={})


class TestReproduceOldBehaviorStorageManager:
    def test_context(self, tmp_path):
        # check that the paths are the ones we expect
        base = tmp_path / "working"
        manager = ReproduceOldBehaviorStorageManager.from_old_args(
            scratch_basedir=base,
            shared_basedir=base
        )
        dag_label = "dag"
        unit_label = "unit"
        expected_scratch = "working/dag/scratch_unit_attempt_0/scratch.txt"
        expected_shared = "working/dag/shared_unit_attempt_0/shared.txt"
        expected_perm = "working/dag/shared_unit_attempt_0/perm.txt"
        with manager.running_dag(dag_label) as dag_ctx:
            with dag_ctx.running_unit(
                dag_label, unit_label, attempt=0
            ) as ctx:
                scratch_f = ctx.scratch / "scratch.txt"
                shared_f = ctx.shared / "shared.txt"
                perm_f = ctx.permanent / "perm.txt"

        found_scratch = pathlib.Path(scratch_f).relative_to(tmp_path)
        found_shared = pathlib.Path(shared_f.fspath).relative_to(tmp_path)
        found_perm = pathlib.Path(perm_f.fspath).relative_to(tmp_path)

        assert str(found_scratch) == expected_scratch
        assert str(found_shared) == expected_shared
        assert str(found_perm) == expected_perm
        # the label is the relative path to the base directory for a
        # FileStorage
        assert "working/" + shared_f.label == expected_shared
        assert "working/" + perm_f.label == expected_perm


@pytest.mark.parametrize('keep_shared', [False, True])
@pytest.mark.parametrize('keep_scratch', [False, True])
def test_execute_dag(tmpdir, keep_shared, keep_scratch, writefile_dag):

    with tmpdir.as_cwd():
        shared = pathlib.Path('shared')
        shared.mkdir(parents=True)

        scratch = pathlib.Path('scratch')
        scratch.mkdir(parents=True)
    
        # run dag
        execute_DAG(writefile_dag,
                    shared_basedir=shared,
                    scratch_basedir=scratch,
                    keep_shared=keep_shared,
                    keep_scratch=keep_scratch)
    
        # check outputs are as expected
        # will have produced 4 files in scratch and shared directory
        for pu in writefile_dag.protocol_units:
            identity = pu.inputs['identity']
            shared_file = (shared
                           / f'shared_{str(pu.key)}_attempt_0'
                           / f'unit_{identity}_shared.txt')
            scratch_file = (scratch
                            / f'scratch_{str(pu.key)}_attempt_0'
                            / f'unit_{identity}_scratch.txt')
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
