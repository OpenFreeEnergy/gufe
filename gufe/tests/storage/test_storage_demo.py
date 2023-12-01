import pytest

import pathlib

import gufe
from gufe.storage.externalresource import MemoryStorage, FileStorage
from gufe.storage.storagemanager import StorageManager
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


class ExecutionStorageDemoTest:
    """
    Template method pattern ABC for tests of StorageDemoProtocol execution.

    Using template method here because it ensures that all aspects get
    tested for all implementations, even though individual aspects may
    differ between different setups.
    """
    def get_shared_and_permanent(self):
        raise NotImplementedError()

    @staticmethod
    def _parse_keep(keep):
        return (
            'scratch' in keep,
            'staging' in keep,
            'shared' in keep,
            'empties' in keep
        )

    def assert_dag_result(self, result, demo_dag, storage_manager):
        """Test that the ProtocolDAGResult has the expected contents.

        This should be preserved across all execution methods.
        """
        u1_label = self.u1_label(demo_dag)
        keep_scratch = storage_manager.keep_scratch

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

    def assert_shared_and_permanent(self, storage_manager, dag):
        """Check the final status of the shared and permanent containers.

        The can depend on the relation between the shared and permanent
        external storage containers. For example, if they are the same
        object, the final contents of permament will also include the final
        contents of shared (and vice versa).

        Default behavior here is for the case of distinct backends.
        """
        shared = storage_manager.shared_root
        permanent = storage_manager.permanent_root
        u1_label = self.u1_label(dag)
        keep_shared = storage_manager.keep_shared

        perm_file = f"{u1_label}/permanent.txt"
        shared_file = f"{u1_label}/shared.txt"

        assert list(permanent.iter_contents()) == [perm_file]
        with permanent.load_stream(perm_file) as f:
            assert f.read() == b"I'm permanent (but I can be shared)"

        if keep_shared:
            assert list(shared.iter_contents()) == [shared_file, perm_file]
            with shared.load_stream(shared_file) as f:
                assert f.read() == b"I can be shared"
            with shared.load_stream(perm_file) as f:
                assert f.read() == b"I'm permanent (but I can be shared)"
        else:
            assert list(shared.iter_contents()) == []

    def assert_scratch(self, storage_manager):
        """Check the final status of the scratch directory.

        This will change if the scratch is within the staging root directory
        (for cases where we want to keep one of staging/scratch and not the
        other; empty directories might get deleted in one case).
        """
        scratch = storage_manager.scratch_root
        keep_scratch = storage_manager.keep_scratch
        del_empty_dirs = storage_manager.delete_empty_dirs
        assert scratch.is_dir()

        if keep_scratch:
            n_expected = 1 if del_empty_dirs else 2
            dag_dir = scratch / "scratch/dag"
            assert len(list(dag_dir.iterdir())) == n_expected
        else:
            assert 'scratch' not in list(scratch.iterdir())

    def assert_staging(self, storage_manager, dag):
        """Check the final status of the staging directory.

        Behavior here will change if staging overlaps with a FileStorage for
        either shared or permanent.
        """
        keep_staging = storage_manager.keep_staging
        u1_label = self.u1_label(dag)
        scratch_root = storage_manager.scratch_root
        u1_staging = scratch_root / ".staging" / u1_label

        if keep_staging:
            assert (u1_staging / "shared.txt").exists()
            assert (u1_staging / "permanent.txt").exists()
        else:
            assert ".staging" not in list(scratch_root.iterdir())

    @staticmethod
    def u1_label(dag):
        """Unit 1 label"""
        return f"dag/{dag.protocol_units[0].key}_attempt_0"

    @staticmethod
    def u2_label(dag):
        """Unit 2 label"""
        return f"dag/{dag.protocol_units[1].key}_attempt_0"

    def get_storage_manager(self, keep, tmp_path):
        keep_scr, keep_sta, keep_sha, empties = self._parse_keep(keep)
        del_empty_dirs = not empties
        shared, permanent = self.get_shared_and_permanent()

        storage_manager = StorageManager(
            scratch_root=tmp_path,
            shared_root=shared,
            permanent_root=permanent,
            keep_scratch=keep_scr,
            keep_staging=keep_sta,
            keep_shared=keep_sha,
            delete_empty_dirs=del_empty_dirs,
        )
        return storage_manager

    @pytest.mark.parametrize('keep', [
        'nothing', 'scratch', 'staging', 'shared', 'scratch,staging',
        'scratch,shared', 'staging,shared', 'scratch,staging,shared',
        'scratch,empties', 'scratch,shared,empties',
    ])
    def test_execute_dag(self, demo_dag, keep, tmp_path):
        storage_manager = self.get_storage_manager(keep, tmp_path)

        dag_label = "dag"
        result = new_execute_DAG(demo_dag, dag_label, storage_manager,
                                 raise_error=True, n_retries=2)

        self.assert_dag_result(result, demo_dag, storage_manager)
        self.assert_shared_and_permanent(storage_manager, demo_dag)
        self.assert_scratch(storage_manager)
        self.assert_staging(storage_manager, demo_dag)


class TestExecuteStorageDemoDiffBackends(ExecutionStorageDemoTest):
    """
    Test execution when permanent and shared are different MemoryStorages.

    This is considered the standard base case; this should be easiest to
    pass, as there should be no special case code that needs to be invoked.
    """
    def get_shared_and_permanent(self):
        return MemoryStorage(), MemoryStorage()


class TestExecuteStorageDemoSameBackend(ExecutionStorageDemoTest):
    """
    Test execution when permanent and shared are the same MemoryStorage.
    """
    def get_shared_and_permanent(self):
        backend = MemoryStorage()
        return backend, backend

    def assert_shared_and_permanent(self, storage_manager, dag):
        shared = storage_manager.shared_root
        permanent = storage_manager.permanent_root
        u1_label = self.u1_label(dag)
        keep_shared = storage_manager.keep_shared

        perm_file = f"{u1_label}/permanent.txt"
        shared_file = f"{u1_label}/shared.txt"

        assert shared is permanent
        # we'll test everything in permanent, because shared is identical

        if keep_shared:
            expected = {perm_file, shared_file}
        else:
            expected = {perm_file}

        assert set(permanent.iter_contents()) == expected
        with permanent.load_stream(perm_file) as f:
            assert f.read() == b"I'm permanent (but I can be shared)"

        if keep_shared:
            with permanent.load_stream(shared_file) as f:
                assert f.read() == b"I can be shared"


class TestExecuteStorageDemoStagingOverlap(TestExecuteStorageDemoSameBackend):
    """
    Test execution when permanent and shared overlap with staging.

    This represents the approach we will probably actually use. In this
    case, we use identical FileStorage for shared and permanent, and those
    overlap with the staging directory. The result is that file locations
    don't actually change.
    """
    def get_shared_and_permanent(self):
        ...  # override the need for this; not the prettiest, but it works

    def get_storage_manager(self, keep, tmp_path):
        keep_scr, keep_sta, keep_sha, empties = self._parse_keep(keep)
        del_empty_dirs = not empties
        backend = FileStorage(tmp_path)
        storage_manager = StorageManager(
            scratch_root=tmp_path,
            shared_root=backend,
            permanent_root=backend,
            keep_scratch=keep_scr,
            keep_staging=keep_sta,
            keep_shared=keep_sha,
            delete_empty_dirs=del_empty_dirs,
            staging="",
        )
        return storage_manager

    def assert_shared_and_permanent(self, storage_manager, dag):
        shared = storage_manager.shared_root
        permanent = storage_manager.permanent_root
        u1_label = self.u1_label(dag)
        keep_shared = storage_manager.keep_shared
        keep_scratch = storage_manager.keep_scratch

        perm_file = f"{u1_label}/permanent.txt"
        shared_file = f"{u1_label}/shared.txt"
        scratch_file = f"scratch/{u1_label}/scratch.txt"

        assert shared is permanent
        # we'll test everything in permanent, because shared is identical

        expected = {perm_file}

        if keep_shared:
            expected.add(shared_file)

        if keep_scratch:
            expected.add(scratch_file)

        assert set(permanent.iter_contents()) == expected
        with permanent.load_stream(perm_file) as f:
            assert f.read() == b"I'm permanent (but I can be shared)"

        if keep_shared:
            with permanent.load_stream(shared_file) as f:
                assert f.read() == b"I can be shared"

        if keep_scratch:
            with permanent.load_stream(scratch_file) as f:
                assert f.read() == b"This is scratch -- can't be shared"

    def assert_staging(self, storage_manager, dag):
        # in this case, keep_staging is ignored in favor of the behavior of
        # keep_shared
        keep_shared = storage_manager.keep_shared
        u1_label = self.u1_label(dag)
        scratch_root = storage_manager.scratch_root
        u1_staging = scratch_root / u1_label

        assert (u1_staging / "permanent.txt").exists()

        if keep_shared:
            assert (u1_staging / "shared.txt").exists()
