import pytest

import os
import pathlib

from gufe.storage.externalresource import MemoryStorage
from gufe.storage.pseudodirectory import SharedRoot


@pytest.fixture
def root(tmp_path):
    external = MemoryStorage()
    external.store_bytes("old_unit/data.txt", b"foo")
    root = SharedRoot(
        scratch=tmp_path,
        external=external,
        prefix="new_unit",
        delete_holding=False
    )
    return root

@pytest.fixture
def root_with_contents(root):
    with open(root / "data.txt", mode='wb') as f:
        f.write(b"bar")

    return root

class TestSharedRoot:
    @pytest.mark.parametrize('pathlist', [
        ['file.txt'], ['dir', 'file.txt']
    ])
    def test_path(self, root, pathlist):
        path = root
        for p in pathlist:
            path = path / p

        inner_path = os.sep.join(pathlist)
        actual_path = root.shared_dir / inner_path

        assert pathlib.Path(path) == actual_path

    def test_read_old(self, root):
        # When the file doesn't exist locally, it should be pulled down the
        # first time that we register the path.

        # initial conditions, without touching SharedRoot/SharedPath
        label = "old_unit/data.txt"
        on_filesystem = root.scratch / root.holding / label
        assert not on_filesystem.exists()
        assert root.external.exists(label)

        # when we create the specific SharedPath, it registers and
        # "downloads" the file
        old_shared = root.get_other_shared_dir("old_unit")
        filepath = old_shared / "data.txt"
        assert pathlib.Path(filepath) == on_filesystem
        assert on_filesystem.exists()

        # let's just be sure we can read in the data as desired
        with open(filepath, mode='rb') as f:
            assert f.read() == b"foo"

    def test_write_new(self, root):
        label = "new_unit/somefile.txt"
        on_filesystem = root.scratch / root.holding / label
        assert not on_filesystem.exists()
        with open(root / "somefile.txt", mode='wb') as f:
            f.write(b"testing")

        # this has been written to disk in scratch, but not yet saved to
        # external storage
        assert on_filesystem.exists()
        assert not root.external.exists(label)

    def test_write_old_fail(self, root):
        old_shared = root.get_other_shared_dir("old_unit")
        with pytest.raises(IOError, match="read-only"):
            old_shared / "foo.txt"

    def test_transfer_to_external(self, root_with_contents):
        path = list(root_with_contents.registry)[0]  # only 1
        assert not root_with_contents.external.exists(path.label)

        root_with_contents.transfer_holding_to_external()
        assert root_with_contents.external.exists(path.label)

        with root_with_contents.external.load_stream(path.label) as f:
            assert f.read() == b"bar"

    def test_transfer_to_external_no_file(self, root):
        ...

    def test_tranfer_to_external_directory(self, root):
        ...

    def test_del(self):
        ...

    def test_existing_local_and_external(self, root):
        ...

    def test_existing_local_and_external_conflict(self, root):
        ...

    def test_no_transfer_for_read_only(self, root):
        ...
