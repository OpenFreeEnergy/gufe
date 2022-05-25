import pytest

from gufe.storage.resultstore import ResultStore

from gufe.storage.externalresource import FileStorage
from gufe.storage.metadatastore import JSONMetadataStore
from gufe.storage.errors import MissingExternalResourceError


@pytest.fixture
def result_store(tmpdir):
    external = FileStorage(tmpdir)
    metadata = JSONMetadataStore(external)
    result_store = ResultStore(external, metadata)
    result_store.store("path/to/foo.txt", "foo".encode('utf-8'))
    return result_store


class TestResultStore:
    def test_store(self, result_store):
        # first check the thing stored during the fixture
        metadata_store = result_store.metadata_store
        foo_loc = "path/to/foo.txt"
        assert len(metadata_store) == 1
        assert foo_loc in metadata_store
        assert result_store.external_store.exists(foo_loc)

        # also explicitly test storing here
        bar_loc = "path/to/bar.txt"
        result_store.store(bar_loc, "bar".encode('utf-8'))
        assert len(metadata_store) == 2
        assert bar_loc in metadata_store
        assert result_store.external_store.exists(bar_loc)
        bar_hash = metadata_store[bar_loc]
        external = result_store.external_store
        with external.load_stream(bar_loc, bar_hash) as f:
            assert f.read().decode('utf-8') == "bar"

    def test_iter(self, result_store):
        assert list(result_store) == ["path/to/foo.txt"]

    def test_find_missing_files(self, result_store):
        result_store.metadata_store.store_metadata("fake/file.txt",
                                                   "1badc0de")

        assert result_store.find_missing_files() == ["fake/file.txt"]

    def test_load_stream(self, result_store):
        with result_store.load_stream('path/to/foo.txt') as f:
            contents = f.read()

        assert contents.decode('utf-8') == "foo"

    def test_load_stream_missing(self, result_store):
        with pytest.raises(MissingExternalResourceError, match="not found"):
            result_store.load_stream("path/does/not/exist.txt")

