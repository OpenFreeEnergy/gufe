import pytest
from gufe.tests.storage.test_storagemanager import dag_units, LifecycleHarness
from gufe.tests.storage.test_storage_demo import demo_dag
from gufe.storage.unitstoragemanager import PerUnitStorageManager
from gufe.tokenization import from_dict

from gufe.storage.externalresource import MemoryStorage, FileStorage
from gufe.protocols.protocoldag import Context, ProtocolDAGResult

import json

# TODO: execute_unit should actually be moved somewhere else; this is likely
# to be the starting point for a real approach to do that
def execute_unit(dag_label, protocolunit, storage_manager, attempt, inputs):
    with storage_manager.running_dag(dag_label) as dag_ctx:
        with dag_ctx.running_unit(
            dag_label,
            protocolunit.key,
            attempt=attempt
        ) as (scratch, shared, perm):
            context = Context(shared=shared,
                              scratch=scratch,
                              permanent=perm)

            unit_result = protocolunit.execute(context=context,
                                               raise_error=False,
                                               **inputs)

    return unit_result

# the next functions will probably become conveniences on the warehouse
def _result_filenames_for_unit(unit, results_dir):
    yield from (
        f for f in results_dir.iterdir()
        if f.name.startswith(unit.key)
    )

def load_results_for_unit(unit, results_dir, storage_manager):
    for filename in _result_filenames_for_unit(unit, results_dir):
        with open(filename, mode='r') as f:
            res = from_dict(json.load(f, cls=storage_manager.json_decoder))

        yield res


def get_inputs(unit, results_dir, storage_manager):
    inputs = {}
    for inp_name, inp_unit in unit.inputs.items():
        for res in load_results_for_unit(inp_unit, results_dir,
                                         storage_manager):
            if res.ok:
                # there should be only 1 unit result that is ok, although
                # we're not being safe about that in this little demo
                inputs[inp_name] = res
                break

    return inputs


def get_attempt_number(unit, results_dir):
    return len(list(_result_filenames_for_unit(unit, results_dir)))


def save_unit_result(result, storage_manager, results_dir, attempt):
    fname = results_dir / f"{result.source_key}_attempt_{attempt}.json"
    dct = result.to_dict()   # real approach should be more efficient
    with open(fname, mode='w') as f:
        f.write(json.dumps(dct, cls=storage_manager.json_encoder))

# now we have a test-only method, which will fake independent processes
# (although the actual executor will have some similar things)
def execute_per_unit(protocoldag, storage_manager, results_dir):
    # fake like we're executing each unit in a different process
    dag_label = "dag"  # is this needed? check with other version
    for num, unit in enumerate(protocoldag.protocol_units):
        # when you run a unit, get its info
        inputs = get_inputs(unit, results_dir, storage_manager)
        attempt = get_attempt_number(unit, results_dir)
        unit_result = execute_unit(dag_label, unit, storage_manager,
                                   attempt, inputs)
        save_unit_result(unit_result, storage_manager, results_dir, attempt)

        # now let's force the unit_result to get cleared from memory
        del unit_result

    # reload all the serialized units (this would be a task generated after
    # the normal tasks to create a DAG result) -- results_dir is specific to
    # this DAG
    unit_results = []
    for fname in results_dir.iterdir():
        with open(fname, mode='r') as f:
            unit_results.append(
                from_dict(json.load(f, cls=storage_manager.json_decoder))
            )

    dag_result = ProtocolDAGResult(
        protocol_units=protocoldag.protocol_units,
        protocol_unit_results=unit_results,
        transformation_key=protocoldag.transformation_key,
        # NOTE: this function doesn't yet allow extends, etc.
    )
    return dag_result


def test_execute_per_unit(tmp_path, demo_dag):
    results = tmp_path / "result_objs"
    scratch = tmp_path / "working"
    shared = MemoryStorage()
    permanent = MemoryStorage()
    results.mkdir(parents=True, exist_ok=True)
    manager = PerUnitStorageManager(
        scratch_root=scratch,
        shared_root=shared,
        permanent_root=permanent,
    )
    dag_result = execute_per_unit(
        demo_dag,
        manager,
        results
    )
    assert dag_result.ok()
    # TODO: further asserts to make sure everything behaved as expected
