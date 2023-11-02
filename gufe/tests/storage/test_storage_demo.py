import pytest

import pathlib

import gufe
from gufe.storage.externalresource import MemoryStorage
from gufe.storage.storagemanager import StorageManager, NewStorageManager
from gufe.storage.stagingdirectory import StagingPath
from gufe.protocols.protocoldag import new_execute_DAG

"""
This module contains complete integration tests for the storage lifecycle,
using an actual protocol as an example.

These tests are largely redundant from the perspective of unit testing, but
the :class:`.StoragedDemoProtocol` is useful as an example for
implementation. Furthermore, as integration tests, they ensure that the
whole setup works together.
"""


class Unit1(gufe.ProtocolUnit):
    def _execute(self, ctx):
        share_file = ctx.shared / "shared.txt"
        with open(share_file, mode='w') as f:
            f.write("I can be shared")

        perm_file = ctx.permanent / "permanent.txt"
        with open(perm_file, mode='w') as f:
            f.write("I'm permanent (but I can be shared)")

        scratch_file = ctx.scratch / "scratch.txt"
        with open(scratch_file, mode='w') as f:
            f.write("This is scratch -- can't be shared")

        return {'share_file': share_file,
                'perm_file': perm_file,
                'scratch_file': scratch_file}


class Unit2(gufe.ProtocolUnit):
    def _execute(self, ctx, unit1_result):
        u1_outputs = unit1_result.outputs

        outputs = {}
        for file_label, file in unit1_result.outputs.items():
            # import pdb; pdb.set_trace()
            # labels are, e.g., share_file; file is StagingPath
            key = f"{file_label}_contents"
            try:
                with open(file, mode='r') as f:
                    outputs[key] = f.read()
            except FileNotFoundError:
                outputs[key] = "File not found"

        return outputs


class StorageDemoProtocol(gufe.Protocol):
    @classmethod
    def _default_settings(cls):
        return {}

    @classmethod
    def _defaults(cls):
        return {}

    def _create(self, stateA, stateB, mapping, extends):
        u1 = Unit1()
        u2 = Unit2(unit1_result=u1)
        return [u1, u2]

    def _gather(self, protocol_dag_results):
        return {}


# TODO: execute_unit should actually be moved somewhere else; this is likely
# to be the starting point for a real approach to do that
def execute_unit(dag_label, protocolunit, storage_manager, inputs):
    label = f"{str(unit.key)}"
    with storage_manager(running_dag(dag_label)) as dag_ctx:
        with dag_ctx.running_unit(label) as (scratch, shared, perm):
            context = Context(shared=shared,
                              scratch=scratch,
                              permanent=perm)

            unit_result = protocolunit.execute(context,
                                               raise_error=False,
                                               **inputs)

    return unit_result


def execute_per_unit(protocoldag, storage_manager, dag_directory):
    # fake like we're executing each unit in a different process
    all_unit_filenames = []
    dag_label = protocoldag.key  # TODO: we can change this
    for num, unit in enumerate(protocoldag.protocol_units):
        unit_result = execute_unit(dag_label, unit, storage_manager)
        fname = dag_directory / f"result_{num}.json"
        # serialize the unit result
        with open(fname, mode='w') as f:
            f.write(json.dumps(unit_result.to_dict(),
                               cls=storage_manager.json_encoder))

        all_unit_filenames.append(fname)

        # now let's force the unit_result to get cleared from memory
        del unit_result
        assert gc.is_finalized(unit_result)

    ... # TODO: make ProtocolDAGResult

@pytest.fixture
def demo_dag(solvated_ligand, solvated_complex):
    transformation = gufe.Transformation(
        solvated_ligand,
        solvated_complex,
        protocol=StorageDemoProtocol(StorageDemoProtocol.default_settings()),
        mapping=None
    )
    dag = transformation.create()
    return dag


def assert_dag_result(result, u1_label, keep_scratch):
    # no matter how you set up storage, the DAGResult and UnitResults from
    # this protocol should be the same
    assert result.ok
    assert len(result.protocol_unit_results) == 2
    res1, res2 = result.protocol_unit_results
    assert set(res1.outputs) == {'share_file', 'perm_file', 'scratch_file'}
    assert isinstance(res1.outputs['scratch_file'], pathlib.Path)
    assert isinstance(res1.outputs['share_file'], StagingPath)
    assert isinstance(res1.outputs['perm_file'], StagingPath)


    if keep_scratch:
        scratch_res2 = "This is scratch -- can't be shared"
    else:
        scratch_res2 = "File not found"

    assert res2.outputs == {
        'share_file_contents': "I can be shared",
        'perm_file_contents': "I'm permanent (but I can be shared)",
        'scratch_file_contents': scratch_res2
    }

def _parse_keep(keep):
    return (
        'scratch' in keep,
        'staging' in keep,
        'shared' in keep
    )



@pytest.mark.parametrize('keep', [
    'nothing', 'scratch', 'staging', 'shared', 'scratch,staging',
    'scratch,shared', 'staging,shared', 'scratch,staging,shared'
])
def test_execute_DAG(demo_dag, tmp_path, keep):
    keep_scratch, keep_staging, keep_shared = _parse_keep(keep)
    dag = demo_dag

    shared = MemoryStorage()
    permanent = MemoryStorage()
    scratch = tmp_path

    storage_manager = NewStorageManager(
        scratch_root=scratch,
        shared_root=shared,
        permanent_root=permanent,
        keep_scratch=keep_scratch,
        keep_shared=keep_shared,
        keep_staging=keep_staging,
    )

    dag_label = "dag"  # currently unused?
    result = new_execute_DAG(dag, dag_label, storage_manager,
                             raise_error=True, n_retries=3)

    # test the ProtocolDAGResult that comes out
    u1_label = f"{dag.protocol_units[0].key}_attempt_0"
    assert_dag_result(result, u1_label, keep_scratch)

    # test the shared/permanent storage
    assert permanent._data == {
        f'{u1_label}/permanent.txt': b"I'm permanent (but I can be shared)",
    }
    if keep_shared:
        assert shared._data == {
            f'{u1_label}/shared.txt': b"I can be shared",
            # f'{u1_label}/permanent.txt': b"I'm permanent (but I can be shared)",
        }
    else:
        assert shared._data == {}

    # test the directories that we generated
    assert scratch.is_dir()

    if keep_scratch:
        assert len(list((scratch / "scratch").iterdir())) == 2
    else:
        assert 'scratch' not in list(scratch.iterdir())

    if keep_staging:
        assert (tmp_path / ".staging" / u1_label / "shared.txt").exists()
        assert (tmp_path / ".staging" / u1_label / "permanent.txt").exists()
    else:
        assert ".staging" not in list(scratch.iterdir())
