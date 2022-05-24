import pytest
import pathlib
import hashlib
from unittest import mock

from gufe.storage.externalresource import FileStorage, MemoryStorage
from gufe.storage.errors import (
    MissingExternalResourceError, ChangedExternalResourceError
)

# input for use with file_storage fixture
BAR_HASH = hashlib.sha256("bar".encode('utf-8')).hexdigest()

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
    def _assert_file_with_contents(self, filename, contents):
        assert pathlib.Path(filename).exists()
        with open(filename, mode='rb') as f:
            assert f.read() == contents.encode('utf-8')

    @pytest.mark.parametrize('filename', [
        'foo.txt', 'notexisting.txt', 'with/directory.txt'
    ])
    def test_exists(self, filename, file_storage):
        expected = {
            'foo.txt': True,
            'notexisting.txt': False,
            'with/directory.txt': True
        }[filename]

        assert file_storage.exists(filename) == expected

    def test_store(self, file_storage):
        fileloc = file_storage.root_dir / "bar.txt"
        assert not fileloc.exists()
        as_bytes = "This is bar".encode('utf-8')
        # mock out the metadata (hash)
        mock_hash = mock.Mock(
            return_value=mock.Mock(
                hexdigest=mock.Mock(return_value="deadbeef")
            )
        )
        with mock.patch('hashlib.sha256', mock_hash):
            loc, sha = file_storage.store("bar.txt", as_bytes)

        assert loc == str(fileloc)
        assert sha == "deadbeef"
        assert fileloc.exists()
        with open(fileloc, 'rb') as f:
            assert as_bytes == f.read()

    def test_delete(self, file_storage):
        path = file_storage.root_dir / "foo.txt"
        assert path.exists()
        file_storage.delete("foo.txt")
        assert not path.exists()

    def test_delete_error_not_existing(self, file_storage):
        with pytest.raises(MissingExternalResourceError,
                           match="does not exist"):
            file_storage.delete("baz.txt")

    def test_get_filename(self, file_storage):
        expected = file_storage.root_dir / "foo.txt"
        filename = file_storage.get_filename("foo.txt", BAR_HASH)
        assert filename == str(expected)

    @pytest.mark.parametrize('as_bytes', [True, False])
    def test_load_stream(self, file_storage, as_bytes):
        with file_storage.load_stream("foo.txt", BAR_HASH, as_bytes) as f:
            results = f.read()
            if as_bytes:
                results = results.decode('utf-8')

            assert results == "bar"

    def test_load_stream_error_missing(self, file_storage):
        with pytest.raises(MissingExternalResourceError):
            file_storage.load_stream("baz.txt", "1badc0de")

    def test_load_stream_error_bad_hash(self, file_storage):
        # for the test, instead of changing the contents, we change the hash
        with pytest.raises(ChangedExternalResourceError):
            file_storage.load_stream("foo.txt", "1badc0de")

    def test_load_stream_allow_bad_hash(self, file_storage):
        allow_changed = file_storage.allow_changed
        file_storage.allow_changed = True
        with pytest.warns(UserWarning, match="Hash mismatch"):
            file = file_storage.load_stream("foo.txt", "1badc0de")

        with file as f:
            assert f.read().decode("utf-8") == "bar"
