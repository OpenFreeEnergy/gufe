import gufe

"""
This module contains a complete integration test for the
"""


class Unit1(gufe.ProtocolUnit):
    def _execute(ctx):
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
    def _execute(ctx, unit1_result):
        u1_outputs = unit1_result.outputs
        share_file = u1_outputs['share_file']
        perm_file = u1_outputs['perm_file']
        scratch_file = u1_outputs['scratch_file']

        outputs = {}
        for file_label, file in unit1_result.outputs.items():
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

    def _create(self, *, stateA, stateB, mapping, extends, name,
                transformation_key):
        u1 = Unit1()
        u2 = Unit2(unit1_results=u1)
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


def test_execute_DAG(solvated_ligand, solvated_complex, tmp_path):
    transformation = gufe.Transformation(
        solvated_liquid,
        solvated_complex,
        protocol=StorageDemoProtocol(),
        mapping=None
    )
    dag = transformation.create()
    shared = MemoryStorage()
    permanent = MemoryStorage()
    storage_manager = StorageManager(
        scratch_root=tmp_path / "scratch",
        shared_root=shared,
        permanent_root=permanent
    )
    result = gufe.protocols.new_execute_DAG(dag, storage_manager,
                                            n_retries=3)
    assert result.ok
    assert len(result.protocol_unit_results) == 2
    res1, res2 = result.protocol_unit_results
    assert set(res1.outputs) == {'share_file', 'perm_file', 'scratch_file'}
    # further tests of res1?

    assert res2.outputs = {
        ...
    }

    assert shared._data = {
        'shared.txt': "I can be shared"
    }
    assert permanent._data = {
        'permanent.txt': "I'm permanent (but I can be shared)",
    }
    # test that scratch is empty




