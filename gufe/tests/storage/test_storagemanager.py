import pytest
from gufe.storage.storagemanager import (
    StorageManager, _storage_path_conflict
)
from gufe.storage.stagingdirectory import StagingDirectory
from gufe.storage.externalresource import MemoryStorage, FileStorage
from pathlib import Path

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

@pytest.fixture
def dag_units():
    class Unit1:
        key = "unit1"

        def run(self, scratch, shared, permanent):
            (scratch / "foo.txt").touch()
            with open(shared / "bar.txt", mode='w') as f:
                f.write("bar was written")
            with open(permanent / "baz.txt", mode='w') as f:
                f.write("baz was written")

            return "done 1"

    class Unit2:
        key = "unit2"

        def run(self, scratch, shared, permanent):
            (scratch / "foo2.txt").touch()
            with shared.other_shared("unit1") as prev_shared:
                with open(prev_shared / "bar.txt", mode='r') as f:
                    bar = f.read()

                # note that you can open a file from permanent as if it was
                # from shared -- everything in permanent is in shared
                with open(prev_shared / "baz.txt", mode='r') as f:
                    baz = f.read()

            return {"bar": bar, "baz": baz}

    return [Unit1(), Unit2()]


@pytest.mark.parametrize('manager', ['std'])
def test_lifecycle(request, manager, dag_units):
    # heavy integration test to ensure that the whole process works
    # this is the primary test of _DAGStorageManager
    storage_manager = request.getfixturevalue(f"storage_manager_{manager}")
    permanent_root = storage_manager.permanent_root
    shared_root = storage_manager.shared_root
    results = []
    unit1_dir = Path(storage_manager.get_shared("dag_label", "unit1"))
    scratch1 = Path(storage_manager.get_scratch("dag_label", "unit1"))
    scratch2 = Path(storage_manager.get_scratch("dag_label", "unit2"))
    barfile = unit1_dir / "bar.txt"
    bazfile = unit1_dir / "baz.txt"
    foofile = scratch1 / "foo.txt"
    foo2file = scratch2 / "foo2.txt"

    all_files = {barfile, bazfile, foofile, foo2file}
    with storage_manager.running_dag("dag_label") as dag_ctx:
        for unit in dag_units:
            with dag_ctx.running_unit(unit) as (scratch, shared, permanent):
                results.append(unit.run(scratch, shared, permanent))

                # check that the expected files are found in staging
                exists = {
                    "unit1": {barfile, bazfile, foofile},
                    "unit2": {foo2file, bazfile}
                }[unit.key]

                for file in exists:
                    assert file.exists()

                for file in all_files - exists:
                    assert not file.exists()

                # check that shared store is as expected
                expected_in_shared = {
                    "unit1": set(),
                    "unit2": {"unit1/bar.txt", "unit1/baz.txt"}
                }[unit.key]
                assert set(shared_root.iter_contents()) == expected_in_shared
                # check that permanent store is empty
                assert list(permanent_root.iter_contents()) == []
            # AFTER THE RUNNING_UNIT CONTEXT
            # Same for both units because unit2 doesn't add anything to
            # shared/permanent
            # Files staged for shared should be transferred to shared and
            # removed from the staging directories; files staged for
            # permanent should remain
            for_permanent = {bazfile}
            for file in for_permanent:
                assert file.exists()

            for file in all_files - for_permanent:
                assert not file.exists()

            # check that we have things in shared
            expected_in_shared = {"unit1/bar.txt", "unit1/baz.txt"}
            assert set(shared_root.iter_contents()) == expected_in_shared
            # ensure that we haven't written to permanent yet
            assert list(permanent_root.iter_contents()) == []
    # AFTER THE RUNNING_DAG CONTEXT
    # all staged files should be deleted
    for file in all_files:
        assert not file.exists()
    # shared still contains everything it had; but this isn't something we
    # guarantee, so we don't actually test for it
    # assert set(shared_root.iter_contents()) == {"unit1/bar.txt",
    #                                             "unit1/baz.txt"}
    assert list(permanent_root.iter_contents()) == ["unit1/baz.txt"]
    # check the results
    assert results == [
        "done 1",
        {"bar": "bar was written", "baz": "baz was written"}
    ]

def test_lifecycle_keep_scratch_and_holding():
    ...

def test_storage_path_conflict_ok(tmp_path):
    # if the filestorage root is not in the given path, no conflict
    external = FileStorage(tmp_path / "foo" / "bar")
    path = tmp_path / "foo" / "baz"
    assert _storage_path_conflict(external, path) is False

def test_storage_path_conflict_not_filestorage(tmp_path):
    # if the external resource isn't a FileStorage, no conflict
    external = MemoryStorage()
    path = tmp_path / "foo" / "baz"
    assert _storage_path_conflict(external, path) is False

def test_storage_path_conflict_problem(tmp_path):
    # if the filestorage root is in the given path, we have a conflict
    external = FileStorage(tmp_path / "foo" / "bar")
    path = tmp_path / "foo"
    assert _storage_path_conflict(external, path) is True


class TestStorageManager:
    def test_get_scratch(self, storage_manager_std):
        scratch = storage_manager_std.get_scratch("dag_label", "unit_label")
        assert str(scratch).endswith("dag_label/scratch/unit_label")
        assert isinstance(scratch, Path)

    def test_get_permanent(self, storage_manager_std):
        perm = storage_manager_std.get_permanent("dag_label", "unit_label")
        assert perm.__fspath__().endswith("dag_label/.holding/unit_label")
        assert isinstance(perm, StagingDirectory)

    def test_get_shared(self, storage_manager_std):
        shared = storage_manager_std.get_shared("dag_label", "unit_label")
        assert shared.__fspath__().endswith("dag_label/.holding/unit_label")
        assert isinstance(shared, StagingDirectory)
