import pytest
from gufe.storage.storagemanager import (
    StorageManager, _storage_path_conflict
)
from gufe.storage.externalresource import MemoryStorage

@pytest.fixture
def storage_manager_std(tmp_path):
    return StorageManager(
        scratch_root=tmp_path / "scratch",
        shared_root=MemoryStorage(),
        permanent_root=MemoryStorage()
    )

@pytest.fixture
def storage_manager_holding_overlaps_shared(tmp_path):
    ...

@pytest.fixture
def storage_manager_holding_overlaps_permanent(tmp_path):
    ...

@pytest.mark.parametrize('manager', ['std'])
def test_lifecycle(request, manager, dag_units):
    # heavy integration test to ensure that the whole process works
    # this is the primary test of _DAGStorageManager
    storage_manager = request.getfixture(f"storage_manager_{manager}")
    with storage_manager.running_dag("dag_label") as dag_ctx:
        for unit in dag_units:
            with dag_ctx.running_unit(unit) as (scratch, shared, permanent):
                results.append(unit.run(scratch, shareed, permanent))
                # TODO: asserts at this level
                # all files exist; are where we expect them
            # TODO: asserts at this level
    # TODO: asserts at this level

def test_lifecycle_keep_scratch_and_holding():
    ...

def test_storage_path_conflict_ok():
    ...

def test_storage_path_conflict_problem():
    ...

class TestStorageManager:
    def test_get_scratch():
        ...

    def test_get_permanent():
        ...

    def test_get_shared():
        ...
