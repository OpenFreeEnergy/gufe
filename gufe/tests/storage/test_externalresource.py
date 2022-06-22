import pytest
import pathlib
import hashlib
import os
from unittest import mock

from gufe.storage.externalresource import FileStorage, MemoryStorage
from gufe.storage.errors import (
    MissingExternalResourceError, ChangedExternalResourceError
)

# NOTE: Tests for the abstract base are just part of the tests of its
# subclasses


@pytest.fixture
def file_storage(tmp_path):
    """Create some files and a FileStorage object.

    Files under tmp_path are:
    * foo.txt : contents "bar"
    * with/directory.txt : contents "in a directory"
    """
    with open(tmp_path / 'foo.txt', 'wb') as foo:
        foo.write("bar".encode("utf-8"))

    inner_dir = tmp_path / 'with'
    inner_dir.mkdir()
    with open(inner_dir / 'directory.txt', 'wb') as with_dir:
        with_dir.write("in a directory".encode("utf-8"))

    return FileStorage(tmp_path)


class TestFileStorage:
    @pytest.mark.parametrize('filename, expected', [
        ('foo.txt', True),
        ('notexisting.txt', False),
        ('with/directory.txt', True),
    ])
    def test_exists(self, filename, expected, file_storage):
        assert file_storage.exists(filename) == expected

    def test_store_bytes(self, file_storage):
        fileloc = file_storage.root_dir / "bar.txt"
        assert not fileloc.exists()
        as_bytes = "This is bar".encode('utf-8')

        file_storage.store_bytes("bar.txt", as_bytes)

        assert fileloc.exists()
        with open(fileloc, 'rb') as f:
            assert as_bytes == f.read()

    def test_store_path(self, file_storage):
        orig_file = file_storage.root_dir / ".hidden" / "bar.txt"
        orig_file.parent.mkdir()
        as_bytes = "This is bar".encode('utf-8')
        with open(orig_file, 'wb') as f:
            f.write(as_bytes)

        fileloc = file_storage.root_dir / "bar.txt"
        assert not fileloc.exists()

        file_storage.store_path(fileloc, orig_file)

        assert fileloc.exists()
        with open(fileloc, 'rb') as f:
            assert as_bytes == f.read()

    def test_delete(self, file_storage):
        path = file_storage.root_dir / "foo.txt"
        assert path.exists()
        file_storage.delete("foo.txt")
        assert not path.exists()

    @pytest.mark.parametrize('prefix,expected', [
        ("", {'foo.txt', 'foo_dir/a.txt', 'foo_dir/b.txt'}),
        ("foo", {'foo.txt', 'foo_dir/a.txt', 'foo_dir/b.txt'}),
        ("foo_dir/", {'foo_dir/a.txt', 'foo_dir/b.txt'}),
        ("foo_dir/a", {'foo_dir/a.txt'}),
        ("foo_dir/a.txt", {'foo_dir/a.txt'}),
        ("baz", set()),
    ])
    def test_iter_contents(self, tmp_path, prefix, expected):
        files = [
            'foo.txt',
            'foo_dir/a.txt',
            'foo_dir/b.txt',
        ]
        for file in files:
            path = tmp_path / file
            path.parent.mkdir(parents=True, exist_ok=True)
            assert not path.exists()
            with open(path, 'wb') as f:
                f.write(b"")

        storage = FileStorage(tmp_path)

        assert set(storage.iter_contents(prefix)) == expected

    def test_delete_error_not_existing(self, file_storage):
        with pytest.raises(MissingExternalResourceError,
                           match="does not exist"):
            file_storage.delete("baz.txt")

    def test_get_filename(self, file_storage):
        expected = file_storage.root_dir / "foo.txt"
        filename = file_storage.get_filename("foo.txt")
        assert filename == str(expected)

    def test_load_stream(self, file_storage):
        with file_storage.load_stream("foo.txt") as f:
            results = f.read().decode('utf-8')

        assert results == "bar"

    def test_load_stream_error_missing(self, file_storage):
        with pytest.raises(MissingExternalResourceError):
            file_storage.load_stream("baz.txt")


class TestMemoryStorage:
    def setup(self):
        self.contents = {'path/to/foo.txt': 'bar'.encode('utf-8')}
        self.storage = MemoryStorage()
        self.storage._data = dict(self.contents)

    @pytest.mark.parametrize('expected', [True, False])
    def test_exists(self, expected):
        path = "path/to/foo.txt" if expected else "path/to/bar.txt"
        assert self.storage.exists(path) is expected

    def test_delete(self):
        # checks internal state
        assert 'path/to/foo.txt' in self.storage._data
        self.storage.delete('path/to/foo.txt')
        assert 'path/to/foo.txt' not in self.storage._data

    def test_delete_error_not_existing(self):
        with pytest.raises(MissingExternalResourceError,
                           match="Unable to delete"):
            self.storage.delete('does/not/exist.txt')

    def test_store_bytes(self):
        storage = MemoryStorage()
        for loc, byte_data in self.contents.items():
            storage.store_bytes(loc, byte_data)

        assert storage._data == self.contents  # internal implementation

    def test_store_path(self, tmp_path):
        storage = MemoryStorage()
        for label, data in self.contents.items():
            path = tmp_path / label
            path.parent.mkdir(parents=True, exist_ok=True)
            with open(path, mode='wb') as f:
                f.write(data)

            storage.store_path(label, path)

        assert storage._data == self.contents

    @pytest.mark.parametrize('prefix,expected', [
        ("", {'foo.txt', 'foo_dir/a.txt', 'foo_dir/b.txt'}),
        ("foo", {'foo.txt', 'foo_dir/a.txt', 'foo_dir/b.txt'}),
        ("foo_dir/", {'foo_dir/a.txt', 'foo_dir/b.txt'}),
        ("foo_dir/a", {'foo_dir/a.txt'}),
        ("foo_dir/a.txt", {'foo_dir/a.txt'}),
        ("baz", set()),
    ])
    def test_iter_contents(self, prefix, expected):
        storage = MemoryStorage()
        storage._data = {
            "foo.txt": b"",
            "foo_dir/a.txt": b"",
            "foo_dir/b.txt": b"",
        }

        assert set(storage.iter_contents(prefix)) == expected


    def test_get_filename(self):
        pytest.skip("Not implemented yet")

    def test_load_stream(self):
        path = "path/to/foo.txt"
        with self.storage.load_stream(path) as f:
            results = f.read().decode('utf-8')

        assert results == "bar"
