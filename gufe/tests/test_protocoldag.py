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
        my_id = inputs['identity']

        with open(os.path.join(ctx.shared, f'unit_{my_id}_shared.txt'), 'w') as out:
            out.write(f'unit {my_id} existed!\n')
        with open(os.path.join(ctx.scratch, f'unit_{my_id}_scratch.txt'), 'w') as out:
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

    def _create(self, stateA,  stateB, mapping, extends=None) -> list[gufe.ProtocolUnit]:
        return [
            WriterUnit(identity=i) for i in range(self.settings.n_repeats)  # type: ignore
        ]

    def _gather(self, results):
        return {}


@pytest.fixture()
def writefile_dag():
    s1 = gufe.ChemicalSystem(components={})
    s2 = gufe.ChemicalSystem(components={})

    p = WriterProtocol(settings=WriterProtocol.default_settings())

    return p.create(stateA=s1, stateB=s2, mapping={})


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
            shared_file = os.path.join(shared,
                                       f'shared_{str(pu.key)}_attempt_0',
                                       f'unit_{identity}_shared.txt')
            scratch_file = os.path.join(scratch,
                                        f'scratch_{str(pu.key)}_attempt_0',
                                        f'unit_{identity}_scratch.txt')
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
