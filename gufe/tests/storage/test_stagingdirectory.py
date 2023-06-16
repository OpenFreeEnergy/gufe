import pytest
from unittest import mock
import logging

import os
import pathlib

from gufe.storage.externalresource import MemoryStorage, FileStorage
from gufe.storage.stagingdirectory import (
    SharedStaging, PermanentStaging, _delete_empty_dirs,
    _safe_to_delete_staging
)

@pytest.fixture
def root(tmp_path):
    external = MemoryStorage()
    external.store_bytes("old_unit/data.txt", b"foo")
    root = SharedStaging(
        scratch=tmp_path,
        external=external,
        prefix="new_unit",
        delete_staging=False
    )
    return root

@pytest.fixture
def root_with_contents(root):
    with open(root / "data.txt", mode='wb') as f:
        f.write(b"bar")

    return root

@pytest.fixture
def read_only_with_overwritten(root_with_contents):
    read_only = SharedStaging(
        scratch=root_with_contents.scratch,
        external=root_with_contents.external,
        prefix="old_unit",
        staging=root_with_contents.staging,
        delete_staging=root_with_contents.delete_staging,
        read_only=True
    )
    filename = pathlib.Path(read_only) / "data.txt"
    assert not filename.exists()
    staged = read_only / "data.txt"
    assert filename.exists()
    with open(staged, mode='w') as f:
        f.write("changed")

    return read_only, staged


def test_safe_to_delete_staging_ok(tmp_path):
    external = FileStorage(tmp_path / "foo")
    prefix = "bar"
    staging = tmp_path / "foo" / "baz"
    assert _safe_to_delete_staging(external, staging, prefix)

def test_safe_to_delete_staging_danger(tmp_path):
    external = FileStorage(tmp_path / "foo")
    prefix = "bar"
    staging = tmp_path / "foo" / "bar" / "baz"
    assert not _safe_to_delete_staging(external, staging, prefix)

def test_safe_to_delete_staging_not_filestorage(tmp_path):
    external = MemoryStorage()
    prefix = "bar"
    staging = tmp_path / "bar"
    assert _safe_to_delete_staging(external, staging, prefix)

def test_delete_empty_dirs(tmp_path):
    base = tmp_path / "tmp"
    paths = [
        base / "foo" / "qux" / "qux.txt",

    ]
    dirs = [
        base / "foo" / "bar" / "baz",
        base / "quux",
    ]
    for directory in dirs:
        directory.mkdir(parents=True, exist_ok=True)

    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()

    _delete_empty_dirs(base)
    for path in paths:
        assert path.exists()

    for directory in dirs:
        assert not directory.exists()

    assert not (base / "foo" / "bar").exists()

@pytest.mark.parametrize('delete_root', [True, False])
def test_delete_empty_dirs_delete_root(tmp_path, delete_root):
    base = tmp_path / "tmp"
    dirs = [
        base / "foo" / "bar" / "baz",
        base / "quux",
    ]
    for directory in dirs:
        directory.mkdir(parents=True, exist_ok=True)

    _delete_empty_dirs(base, delete_root=delete_root)

    for directory in dirs:
        assert not directory.exists()

    assert not (base / "foo" / "bar").exists()
    assert base.exists() is not delete_root




class TestSharedStaging:
    def test_repr(self, root):
        r = repr(root)
        assert r.startswith("SharedStaging")
        assert "MemoryStorage" in r
        assert r.endswith(", new_unit)")

    @pytest.mark.parametrize('pathlist', [
        ['file.txt'], ['dir', 'file.txt']
    ])
    def test_path(self, root, pathlist):
        path = root
        for p in pathlist:
            path = path / p

        inner_path = os.sep.join(pathlist)
        actual_path = root.staging_dir / inner_path

        assert pathlib.Path(path) == actual_path

    def test_read_old(self, root):
        # When the file doesn't exist locally, it should be pulled down the
        # first time that we register the path.

        # initial conditions, without touching StagingDirectory/StagingPath
        label = "old_unit/data.txt"
        on_filesystem = root.scratch / root.staging / "old_unit/data.txt"
        assert not on_filesystem.exists()
        assert root.external.exists(label)

        # when we create the specific StagingPath, it registers and
        # "downloads" the file
        old_staging = root.get_other_shared("old_unit")
        filepath = old_staging / "data.txt"
        assert pathlib.Path(filepath) == on_filesystem
        assert on_filesystem.exists()

        # let's just be sure we can read in the data as desired
        with open(filepath, mode='rb') as f:
            assert f.read() == b"foo"

    def test_write_new(self, root):
        label = "new_unit/somefile.txt"
        on_filesystem = root.scratch / root.staging / "new_unit/somefile.txt"
        assert not on_filesystem.exists()
        with open(root / "somefile.txt", mode='wb') as f:
            f.write(b"testing")

        # this has been written to disk in scratch, but not yet saved to
        # external storage
        assert on_filesystem.exists()
        assert not root.external.exists(label)

    def test_write_old_fail(self, root):
        old_staging = root.get_other_shared("old_unit")
        with pytest.raises(IOError, match="read-only"):
            old_staging / "foo.txt"

    def test_transfer_to_external(self, root_with_contents):
        path = list(root_with_contents.registry)[0]  # only 1
        assert not root_with_contents.external.exists(path.label)

        root_with_contents.transfer_staging_to_external()
        assert root_with_contents.external.exists(path.label)

        with root_with_contents.external.load_stream(path.label) as f:
            assert f.read() == b"bar"

    def test_transfer_to_external_no_file(self, root, caplog):
        with mock.patch.object(root, 'register_path'):
            nonfile = root / "does_not_exist.txt"
        # ensure that we've set this up correctly
        assert nonfile not in root.registry
        logger_name = "gufe.storage.stagingdirectory"
        caplog.set_level(logging.INFO, logger=logger_name)
        root.transfer_single_file_to_external(nonfile)
        assert len(caplog.records) == 1
        record = caplog.records[0]
        assert "nonexistent" in record.msg

    def test_transfer_to_external_directory(self, root, caplog):
        directory = root / "directory"
        with open(directory / "file.txt", mode='w') as f:
            f.write("foo")

        logger_name = "gufe.storage.stagingdirectory"
        caplog.set_level(logging.DEBUG, logger=logger_name)
        root.transfer_single_file_to_external(directory)
        assert len(caplog.records) == 1
        record = caplog.records[0]
        assert "Found directory" in record.msg
        assert "not transfering" in record.msg

    def test_single_file_transfer_read_only(self,
                                            read_only_with_overwritten,
                                            caplog):
        read_only, staged = read_only_with_overwritten
        with read_only.external.load_stream("old_unit/data.txt") as f:
            old_contents = f.read()

        assert old_contents == b"foo"
        logger_name = "gufe.storage.stagingdirectory"
        caplog.set_level(logging.DEBUG, logger=logger_name)
        read_only.transfer_single_file_to_external(staged)
        assert len(caplog.records) == 1
        record = caplog.records[0]
        assert "Read-only:" in record.msg
        with read_only.external.load_stream("old_unit/data.txt") as f:
            new_contents = f.read()
        assert old_contents == new_contents

    def test_transfer_read_only(self, read_only_with_overwritten, caplog):
        read_only, staged = read_only_with_overwritten
        with read_only.external.load_stream("old_unit/data.txt") as f:
            old_contents = f.read()

        assert old_contents == b"foo"
        logger_name = "gufe.storage.stagingdirectory"
        caplog.set_level(logging.DEBUG, logger=logger_name)
        read_only.transfer_staging_to_external()
        assert len(caplog.records) == 1
        record = caplog.records[0]
        assert "Read-only:" in record.msg
        with read_only.external.load_stream("old_unit/data.txt") as f:
            new_contents = f.read()
        assert old_contents == new_contents

    def test_cleanup(self, root_with_contents):
        root_with_contents.delete_staging = True  # slightly naughty
        path = pathlib.Path(root_with_contents.__fspath__()) / "data.txt"
        assert path.exists()
        root_with_contents.cleanup()
        assert not path.exists()

    def test_register_cleanup_preexisting_file(self, root):
        filename = pathlib.Path(root.__fspath__()) / "foo.txt"
        filename.touch()
        root.external.store_bytes("new_unit/foo.txt", b"")
        assert len(root.registry) == 0
        assert len(root.preexisting) == 0
        staging = root / "foo.txt"
        assert staging.label == "new_unit/foo.txt"
        assert len(root.registry) == 1
        assert len(root.preexisting) == 1

        assert filename.exists()
        root.cleanup()
        assert filename.exists()


class TestPermanentStage:
    def test_delete_staging_safe(self):
        ...
