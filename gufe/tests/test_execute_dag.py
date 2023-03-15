# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import os
import pathlib
import pytest
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


class WriterSettings(gufe.settings.ProtocolSettings):
    n_repeats: int


class WriterProtocol(gufe.Protocol):
    result_cls = None

    @classmethod
    def _default_settings(cls) -> gufe.Settings:
        return gufe.Settings(
            settings_version=1,
            thermo_settings=gufe.settings.ThermoSettings(),
            forcefield_settings=gufe.settings.OpenMMSystemGeneratorFFSettings(),
            protocol_settings=WriterSettings(n_repeats=4),
        )

    @classmethod
    def _defaults(cls):
        return {}

    def _create(self, stateA,  stateB, mapping=None, extends=None) -> list[gufe.ProtocolUnit]:
        return [
            WriterUnit(identity=i) for i in range(self.settings.protocol_settings.n_repeats)
        ]

    def _gather(self, results):
        return {}


@pytest.fixture()
def writefile_dag():
    s1 = gufe.ChemicalSystem(components={})
    s2 = gufe.ChemicalSystem(components={})

    p = WriterProtocol(settings=WriterProtocol.default_settings())

    return p.create(stateA=s1, stateB=s2, mapping={})


@pytest.fixture
def in_tmpdir(tmpdir):
    with tmpdir.as_cwd():
        yield str(tmpdir)


@pytest.mark.parametrize('shared', [None, 'myshared'])
@pytest.mark.parametrize('scratch', [None, 'myscratch'])
@pytest.mark.parametrize('use_pathlib', [False, True])
@pytest.mark.parametrize('keep_scratch', [False, True])
def test_execute_dag(in_tmpdir, shared, scratch, use_pathlib, keep_scratch, writefile_dag):
    if scratch is None and keep_scratch:
        pytest.skip("Nonsensical combination")

    if use_pathlib:
        if shared is not None:
            shared = pathlib.Path(shared)
        if scratch is not None:
            scratch = pathlib.Path(scratch)
    if shared is not None:
        os.makedirs(shared)
        shared_ = shared
    else:
        shared_ = in_tmpdir
    if scratch is not None:
        os.makedirs(scratch)
        scratch_ = scratch

    # run dag
    execute_DAG(writefile_dag,
                shared=shared,
                scratch=scratch,
                keep_scratch=keep_scratch)

    # check outputs are as expected
    # will have produced 3 files in scratch and shared directory
    assert os.path.exists(os.path.join(shared_, 'unit_0_shared.txt'))
    if scratch is not None:
        if keep_scratch:
            assert os.path.exists(os.path.join(scratch_, 'unit_0_scratch.txt'))
        else:
            assert not os.path.exists(os.path.join(scratch_, 'unit_0_scratch.txt'))

