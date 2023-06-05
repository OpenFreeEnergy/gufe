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
        scratch_root=tmp_path / "working",
        shared_root=MemoryStorage(),
        permanent_root=MemoryStorage()
    )

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
            # TODO: this will change; the inputs should include a way to get
            # the previous shared unit label
            with shared.other_shared("dag/unit1") as prev_shared:
                with open(prev_shared / "bar.txt", mode='r') as f:
                    bar = f.read()

                # note that you can open a file from permanent as if it was
                # from shared -- everything in permanent is in shared
                with open(prev_shared / "baz.txt", mode='r') as f:
                    baz = f.read()

            return {"bar": bar, "baz": baz}

    return [Unit1(), Unit2()]

class LifecycleHarness:
    @pytest.fixture
    def storage_manager(self, tmp_path):
        raise NotImplementedError()

    @staticmethod
    def get_files_dict(storage_manager):
        root = storage_manager.scratch_root
        holding = storage_manager.holding
        return {
            "foo": root / "dag/unit1/scratch/foo.txt",
            "foo2": root / "dag/unit2/scratch/foo2.txt",
            "bar": root / "dag/unit1" / holding / "bar.txt",
            "baz": root / "dag/unit1" / holding / "baz.txt",
        }

    def test_lifecycle(self, storage_manager, dag_units, tmp_path):
        results = []
        dag_label = "dag"
        with storage_manager.running_dag(dag_label) as dag_ctx:
            for unit in dag_units:
                label = f"{dag_ctx.dag_label}/{unit.key}"
                with dag_ctx.running_unit(label) as (scratch, shared, perm):
                    results.append(unit.run(scratch, shared, perm))
                    self.in_unit_asserts(storage_manager, label)
                self.after_unit_asserts(storage_manager, label)
        self.after_dag_asserts(storage_manager)
        assert results == [
            "done 1",
            {"bar": "bar was written", "baz": "baz was written"}
        ]

    def _in_unit_existing_files(self, unit_label):
        raise NotImplementedError()

    def _after_unit_existing_files(self, unit_label):
        raise NotImplementedError()

    def _after_dag_existing_files(self):
        raise NotImplementedError()

    @staticmethod
    def assert_existing_files(files_dict, existing):
        for file in existing:
            assert files_dict[file].exists()

        for file in set(files_dict) - existing:
            assert not files_dict[file].exists()

    def in_unit_asserts(self, storage_manager, unit_label):
        # check that shared and permanent are correct
        shared_root = storage_manager.shared_root
        permanent_root = storage_manager.permanent_root
        expected_in_shared = {
            "dag/unit1": set(),
            "dag/unit2": {"dag/unit1/bar.txt", "dag/unit1/baz.txt"}
        }[unit_label]
        assert set(shared_root.iter_contents()) == expected_in_shared

        assert list(permanent_root.iter_contents()) == []

        # manager-specific check for files
        files_dict = self.get_files_dict(storage_manager)
        existing = self._in_unit_existing_files(unit_label)
        self.assert_existing_files(files_dict, existing)

    def after_unit_asserts(self, storage_manager, unit_label):
        shared_root = storage_manager.shared_root
        permanent_root = storage_manager.permanent_root
        # these are independent of unit label
        expected_in_shared = {"dag/unit1/bar.txt", "dag/unit1/baz.txt"}
        assert set(shared_root.iter_contents()) == expected_in_shared
        assert list(permanent_root.iter_contents()) == []

        files_dict = self.get_files_dict(storage_manager)
        existing = self._after_unit_existing_files(unit_label)
        self.assert_existing_files(files_dict, existing)

    def after_dag_asserts(self, storage_manager):
        permanent_root = storage_manager.permanent_root
        # shared still contains everything it had; but this isn't something
        # we guarantee, so we don't actually test for it:
        # shared_root = storage_manager.shared_root
        # assert set(shared_root.iter_contents()) == {"unit1/bar.txt",
        #                                             "unit1/baz.txt"}
        assert list(permanent_root.iter_contents()) == ["dag/unit1/baz.txt"]

        files_dict = self.get_files_dict(storage_manager)
        existing = self._after_dag_existing_files()
        self.assert_existing_files(files_dict, existing)


class TestStandardStorageManager(LifecycleHarness):
    @pytest.fixture
    def storage_manager(self, storage_manager_std):
        return storage_manager_std

    def _in_unit_existing_files(self, unit_label):
        return {
            "dag/unit1": {'bar', 'baz', 'foo'},
            "dag/unit2": {'foo2', 'baz'}
        }[unit_label]

    def _after_unit_existing_files(self, unit_label):
        # Same for both units because unit2 doesn't add anything to
        # shared/permanent; in this one, only files staged for permanent
        # should remain
        return {'baz'}

    def _after_dag_existing_files(self):
        return set()


def test_lifecycle_keep_scratch_and_holding():
    ...

def test_lifecycle_holding_overlaps_shared(tmp_path):
    ...

def test_lifecycle_holding_overlaps_permanent(tmp_path):
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
        scratch = storage_manager_std.get_scratch("dag_label/unit_label")
        assert str(scratch).endswith("dag_label/unit_label/scratch")
        assert isinstance(scratch, Path)

    def test_get_permanent(self, storage_manager_std):
        perm = storage_manager_std.get_permanent("dag_label/unit_label")
        assert perm.__fspath__().endswith("dag_label/unit_label/.holding")
        assert isinstance(perm, StagingDirectory)

    def test_get_shared(self, storage_manager_std):
        shared = storage_manager_std.get_shared("dag_label/unit_label")
        assert shared.__fspath__().endswith("dag_label/unit_label/.holding")
        assert isinstance(shared, StagingDirectory)
