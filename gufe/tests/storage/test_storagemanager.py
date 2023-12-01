import pytest
from gufe.storage.storagemanager import StorageManager
from gufe.storage.stagingdirectory import StagingDirectory
from gufe.storage.externalresource import MemoryStorage, FileStorage
from pathlib import Path


@pytest.fixture
def storage_manager_std(tmp_path):
    return StorageManager(
        scratch_root=tmp_path / "working",
        shared_root=MemoryStorage(),
        permanent_root=MemoryStorage(),
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
            with (
                shared.root.other_shared("dag/unit1_attempt_0") as prev_shared
            ):
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
        staging = storage_manager.staging
        return {
            "foo": root / "scratch/dag/unit1_attempt_0/foo.txt",
            "foo2": root / "scratch/dag/unit2_attempt_0/foo2.txt",
            "bar": root / staging / "dag/unit1_attempt_0/bar.txt",
            "baz": root / staging / "dag/unit1_attempt_0/baz.txt",
        }

    def test_lifecycle(self, storage_manager, dag_units, tmp_path):
        results = []
        dag_label = "dag"
        with storage_manager.running_dag(dag_label) as dag_ctx:
            for unit in dag_units:
                label = f"{dag_label}/{unit.key}"
                with dag_ctx.running_unit(dag_label, unit.key, attempt=0) as (
                    scratch, shared, perm
                ):
                    results.append(unit.run(scratch, shared, perm))
                    # import pdb; pdb.set_trace()
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

    def _in_staging_shared(self, unit_label, in_after):
        """
        This is to include things when a shared staging directory reports
        that files exist in it.
        """
        return set()

    def _in_staging_permanent(self, unit_label, in_after):
        """
        This is to include things when a permanent staging directory reports
        that files exist in it.
        """
        return set()

    def in_unit_asserts(self, storage_manager, unit_label):
        # check that shared and permanent are correct
        shared_root = storage_manager.shared_root
        permanent_root = storage_manager.permanent_root
        expected_in_shared = {
            "dag/unit1": set(),
            "dag/unit2": {"dag/unit1_attempt_0/bar.txt",
                          "dag/unit1_attempt_0/baz.txt"}
        }[unit_label] | self._in_staging_shared(unit_label, "in")
        assert set(shared_root.iter_contents()) == expected_in_shared

        expected_in_permanent = self._in_staging_permanent(unit_label, "in")
        assert set(permanent_root.iter_contents()) == expected_in_permanent

        # manager-specific check for files
        files_dict = self.get_files_dict(storage_manager)
        existing = self._in_unit_existing_files(unit_label)
        self.assert_existing_files(files_dict, existing)

    def after_unit_asserts(self, storage_manager, unit_label):
        shared_root = storage_manager.shared_root
        permanent_root = storage_manager.permanent_root
        shared_extras = self._in_staging_shared(unit_label, "after")
        permanent_extras = self._in_staging_permanent(unit_label, "after")
        expected_in_shared = {"dag/unit1_attempt_0/bar.txt",
                              "dag/unit1_attempt_0/baz.txt"}
        expected_in_shared |= shared_extras
        assert set(shared_root.iter_contents()) == expected_in_shared
        assert set(permanent_root.iter_contents()) == permanent_extras

        # manager-specific check for files
        files_dict = self.get_files_dict(storage_manager)
        existing = self._after_unit_existing_files(unit_label)
        self.assert_existing_files(files_dict, existing)

    def after_dag_asserts(self, storage_manager):
        permanent_root = storage_manager.permanent_root
        permanent_extras = self._in_staging_permanent('dag/unit2', "after")
        # shared still contains everything it had; but this isn't something
        # we guarantee, so we don't actually test for it, but we could with
        # this:
        # shared_root = storage_manager.shared_root
        # shared_extras = self._in_staging_shared('dag/unit2', "after")
        # expected_in_shared = {"dag/unit1/bar.txt", "dag/unit1/baz.txt"}
        # expected_in_shared |= shared_extras
        # assert set(shared_root.iter_contents()) == expected_in_shared
        expected_in_permanent = ({"dag/unit1_attempt_0/baz.txt"}
                                 | permanent_extras)
        assert set(permanent_root.iter_contents()) == expected_in_permanent

        # manager-specific check for files
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


class TestKeepScratchAndStagingStorageManager(LifecycleHarness):
    @pytest.fixture
    def storage_manager(self, tmp_path):
        return StorageManager(
            scratch_root=tmp_path / "working",
            shared_root=MemoryStorage(),
            permanent_root=MemoryStorage(),
            keep_scratch=True,
            keep_staging=True
        )

    @staticmethod
    def files_after_unit(unit_label):
        unit1 = {'bar', 'baz', 'foo'}
        unit2 = {'foo2', 'baz'}
        return {
            'dag/unit1': unit1,
            'dag/unit2': unit1 | unit2
        }[unit_label]

    def _in_unit_existing_files(self, unit_label):
        return self.files_after_unit(unit_label)

    def _after_unit_existing_files(self, unit_label):
        return self.files_after_unit(unit_label)

    def _after_dag_existing_files(self):
        return self.files_after_unit('dag/unit2')


class TestStagingOverlapsSharedStorageManager(LifecycleHarness):
    @pytest.fixture
    def storage_manager(self, tmp_path):
        root = tmp_path / "working"
        return StorageManager(
            scratch_root=root,
            shared_root=FileStorage(root),
            permanent_root=MemoryStorage(),
            staging="",
        )

    def _in_unit_existing_files(self, unit_label):
        return {
            "dag/unit1": {'foo', 'bar', 'baz'},
            "dag/unit2": {'foo2', 'bar', 'baz'},
        }[unit_label]

    def _after_unit_existing_files(self, unit_label):
        # same for both; all files come from unit 1
        return {"bar", "baz"}

    def _after_dag_existing_files(self):
        # these get deleted because we don't keep shared here
        return set()

    def _in_staging_shared(self, unit_label, in_after):
        bar = "dag/unit1_attempt_0/bar.txt"
        baz = "dag/unit1_attempt_0/baz.txt"
        foo = "scratch/dag/unit1_attempt_0/foo.txt"
        foo2 = "scratch/dag/unit2_attempt_0/foo2.txt"
        return {
            ("dag/unit1", "in"): {bar, baz, foo},
            ("dag/unit1", "after"): {bar, baz},
            ("dag/unit2", "in"): {bar, baz, foo2},
            ("dag/unit2", "after"): {baz}
        }[unit_label, in_after]


class TestStagingOverlapsPermanentStorageManager(LifecycleHarness):
    @pytest.fixture
    def storage_manager(self, tmp_path):
        root = tmp_path / "working"
        return StorageManager(
            scratch_root=root,
            permanent_root=FileStorage(root),
            shared_root=MemoryStorage(),
            staging="",
        )

    def _in_unit_existing_files(self, unit_label):
        return {
            "dag/unit1": {'foo', 'bar', 'baz'},
            "dag/unit2": {"foo2", "baz"},  # no bar because it was temporary
        }[unit_label]

    def _after_dag_existing_files(self):
        return {"baz"}

    def _in_staging_permanent(self, unit_label, in_after):
        bar = "dag/unit1_attempt_0/bar.txt"
        baz = "dag/unit1_attempt_0/baz.txt"
        foo = "scratch/dag/unit1_attempt_0/foo.txt"
        foo2 = "scratch/dag/unit2_attempt_0/foo2.txt"
        return {
            ("dag/unit1", "in"): {bar, baz, foo},
            ("dag/unit1", "after"): {baz},
            ("dag/unit2", "in"): {baz, foo2},
            ("dag/unit2", "after"): {baz}
        }[unit_label, in_after]

    def _after_unit_existing_files(self, unit_label):
        # same for both; all files come from unit 1
        return {"baz"}
