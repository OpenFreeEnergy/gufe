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
        assert r.startswith("StagingDirectory")
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

    @mock.patch.object(SharedStaging, 'register_path')
    def test_transfer_to_external_no_file(self, root, caplog):
        nonfile = root / "does_not_exist.txt"
        # ensure that we've set this up correctly
        assert nonfile not in root.registry
        caplog.set_level(logging.INFO, logger="gufe.storage")
        root.transfer_single_file_to_external(nonfile)
        assert len(caplog.records) == 1




        ...

    def test_tranfer_to_external_directory(self, root):
        ...

    def test_existing_local_and_external(self, root):
        ...

    def test_existing_local_and_external_conflict(self, root):
        ...

    def test_no_transfer_for_read_only(self, root):
        ...
