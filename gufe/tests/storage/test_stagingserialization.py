import pytest
from gufe.storage.stagingserialization import StagingPathSerialization

from gufe.storage.stagingregistry import StagingPath
from gufe.storage.storagemanager import StorageManager
from gufe.storage.externalresource import MemoryStorage, FileStorage

import json
import pathlib
import shutil


@pytest.fixture
def storage_manager(tmp_path):
    return StorageManager(
        scratch_root=tmp_path / "working",
        shared_root=MemoryStorage(),
        permanent_root=MemoryStorage(),
    )


@pytest.fixture
def shared_path(storage_manager):
    label = storage_manager.make_label("dag", "unit", attempt=0)
    path = storage_manager.shared_staging / label / "file.txt"
    with open(path, mode='w') as f:
        f.write("contents here")

    storage_manager.shared_staging.transfer_staging_to_external()
    return path


@pytest.fixture
def permanent_path(storage_manager):
    label = storage_manager.make_label("dag", "unit", attempt=0)
    path = storage_manager.permanent_staging / label / "file.txt"
    with open(path, mode='w') as f:
        f.write("contents here")

    storage_manager.permanent_staging.transfer_staging_to_external()
    return path


@pytest.fixture
def scratch_path(storage_manager):
    scratch_dir = storage_manager._scratch_loc("dag", "unit", attempt=0)
    path = scratch_dir / "file.txt"
    return path


@pytest.fixture
def serialization_handler(storage_manager):
    return StagingPathSerialization(storage_manager)


class TestStagingPathSerialization:
    @pytest.mark.parametrize('pathtype', ['scratch', 'shared', 'permanent'])
    def test_round_trip(self, serialization_handler, pathtype, request):
        # NB: scratch is a pathlib.Path, not a StagingPath. It is tested
        # here to ensure round-trips as part of the overall user story for
        # this, but it doesn't invoke the machinery of the
        # StagingPathSerialization object
        path = request.getfixturevalue(f"{pathtype}_path")
        as_json = json.dumps(
            path,
            cls=serialization_handler.json_handler.encoder
        )
        reloaded = json.loads(
            as_json,
            cls=serialization_handler.json_handler.decoder
        )

        assert path == reloaded

    @pytest.mark.parametrize('pathtype', ['shared', 'permanent'])
    def test_to_dict(self, serialization_handler, pathtype, request):
        path = request.getfixturevalue(f"{pathtype}_path")
        dct = serialization_handler.to_dict(path)
        assert dct == {
            ':container:': pathtype,
            ':label:': "dag/unit_attempt_0/file.txt",
        }

    # tests for specific user stories
    @pytest.mark.parametrize('pathtype', ['shared', 'permanent'])
    def test_reload_file_contents(self, pathtype, request):
        # USER STORY: I am loading my results object, and I will want to use
        # the associated files. This should be transparent, regardless of
        # where the storage is located (e.g., not local storage). (This is
        # actually a test of the staging tools, but we include it here for
        # completeness of the user stories.)
        path = request.getfixturevalue(f"{pathtype}_path")

        # remove the file (remains in the MemoryStorage)
        p = path.as_path()
        assert p.exists()
        p.unlink()
        assert not p.exists()

        # reload the file (NB: nothing special done here; download is
        # transparent to user)
        with open(path, mode='r') as f:
            contents = f.read()

        assert p.exists()
        assert contents == "contents here"

    @pytest.mark.parametrize('pathtype', ['shared', 'permanent'])
    def test_load_results_object_file_not_downloaded(self,
                                                     serialization_handler,
                                                     pathtype, request):
        # USER STORY: I am loading my results object, but I do not need the
        # large stored files. I do not want to download them when they
        # aren't needed.
        path = request.getfixturevalue(f"{pathtype}_path")
        # serialize the path object
        json_str = json.dumps(path, cls=serialization_handler.encoder)

        # delete the path from the directory
        p = pathlib.Path(path)
        assert p.exists()
        p.unlink()
        assert not p.exists()

        # reload the serialized form of the object
        reloaded = json.loads(json_str, cls=serialization_handler.decoder)

        # check that the deserialized version has the path, but that the
        # path does not exist on the filesystem
        assert isinstance(reloaded, StagingPath)
        assert reloaded.label == path.label
        assert reloaded.path == path.path
        assert not p.exists()
        # NOTE: as soon as you call `__fspath__`, the file will download

    @pytest.mark.parametrize('move', ['relative', 'absolute'])
    def test_permanent_storage_moved(self, move, tmp_path, monkeypatch):
        # USER STORY: My permanent storage was a directory on my file
        # system, but I have moved that directory (with use cases of (a) I
        # moved the absolute path; (b) it is at a different relative path
        # with respect to my pwd).
        monkeypatch.chdir(tmp_path)
        old_manager = StorageManager(
            scratch_root="old/scratch",
            shared_root=FileStorage("old/shared"),
            permanent_root=FileStorage("old/permanent")
        )
        old_handler = StagingPathSerialization(old_manager)
        old_path = old_manager.permanent_staging / "dag/unit/result.txt"
        with open(old_path, mode='w') as f:
            f.write("contents here")

        old_manager.permanent_staging.transfer_staging_to_external()
        perm_p = pathlib.Path(tmp_path / "old/permanent/dag/unit/result.txt")
        assert perm_p.exists()

        # serialize the path object
        json_str = json.dumps(old_path, cls=old_handler.encoder)

        # move the storage subdirectory; create a new, associated storage
        # manager/serialization handler
        if move == "relative":
            # change to within t
            monkeypatch.chdir(tmp_path / "old")
            new_manager = StorageManager(
                scratch_root="scratch",
                shared_root=FileStorage("shared"),
                permanent_root=FileStorage("permanent")
            )
            expected_path = tmp_path / "old/permanent/dag/unit/result.txt"
        elif move == "absolute":
            shutil.move(tmp_path / "old", tmp_path / "new")
            new_manager = StorageManager(
                scratch_root="new/scratch",
                shared_root=FileStorage("new/shared"),
                permanent_root=FileStorage("new/permanent")
            )
            expected_path = tmp_path / "new/permanent/dag/unit/result.txt"
        else:  # -no-cov-
            raise RuntimeWarning(f"Bad test parameter '{move}': should be "
                                 "'relative' or 'absolute'")

        new_handler = StagingPathSerialization(new_manager)

        # deserialize the path using the new serialization handler
        reloaded = json.loads(json_str, cls=new_handler.decoder)

        # ensure that the path exists and that the data can be reloaded
        assert isinstance(reloaded, StagingPath)
        assert reloaded.label == old_path.label
        assert pathlib.Path(expected_path).exists()

        with open(reloaded, mode='r') as f:
            contents = f.read()

        assert contents == "contents here"

    def test_two_different_permanent_storages(self, tmp_path):
        # I'm working with files from two different permanent storages. I
        # need to be able to load from both in the same Python process.
        # (NOTE: this user story is primarily to prevent us from changing to
        # a solution based on global/class vars to set context.)
        manager1 = StorageManager(
            scratch_root=tmp_path / "working1",
            shared_root=MemoryStorage(),
            permanent_root=MemoryStorage(),
        )
        manager2 = StorageManager(
            scratch_root=tmp_path / "working2",
            shared_root=MemoryStorage(),
            permanent_root=MemoryStorage(),
        )
        handler1 = StagingPathSerialization(manager1)
        handler2 = StagingPathSerialization(manager2)

        path1 = manager1.permanent_staging / "file1.txt"
        with open(path1, mode='w') as f:
            f.write("contents 1")
        manager1.permanent_staging.transfer_staging_to_external()

        path2 = manager2.permanent_staging / "file2.txt"
        with open(path2, mode='w') as f:
            f.write("contents 2")
        manager2.permanent_staging.transfer_staging_to_external()

        # serialize the paths
        json_str1 = json.dumps(path1, cls=handler1.encoder)
        json_str2 = json.dumps(path2, cls=handler2.encoder)

        # delete all staged files
        assert path1.as_path().exists()
        manager1.permanent_staging.cleanup()
        assert not path1.as_path().exists()

        assert path2.as_path().exists()
        manager2.permanent_staging.cleanup()
        assert not path2.as_path().exists()

        # reload and check contents of both permanent files
        reloaded1 = json.loads(json_str1, cls=handler1.decoder)
        reloaded2 = json.loads(json_str2, cls=handler2.decoder)

        assert isinstance(reloaded1, StagingPath)
        assert reloaded1.label == path1.label
        assert not reloaded1.as_path().exists()
        with open(reloaded1, mode='r') as f:
            assert f.read() == "contents 1"

        assert isinstance(reloaded2, StagingPath)
        assert reloaded2.label == path2.label
        assert not reloaded2.as_path().exists()
        with open(reloaded2, mode='r') as f:
            assert f.read() == "contents 2"
