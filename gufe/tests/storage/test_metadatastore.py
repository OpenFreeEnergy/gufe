import pytest
import json

from gufe.storage.metadatastore import JSONMetadataStore
from gufe.storage.externalresource import FileStorage
from gufe.storage.externalresource.base import Metadata
from gufe.storage.errors import MissingExternalResourceError


@pytest.fixture
def json_metadata(tmpdir):
    metadata_dict = {'path/to/foo.txt': {'md5': 'bar'}}
    external_store = FileStorage(str(tmpdir))
    with open(tmpdir / 'metadata.json', mode='wb') as f:
        f.write(json.dumps(metadata_dict).encode('utf-8'))
    json_metadata = JSONMetadataStore(external_store)
    return json_metadata


class TestJSONMetadataStore:
    def test_store_metadata(self, json_metadata):
        meta = Metadata(md5="other")
        json_metadata.store_metadata("path/to/other.txt", meta)
        base_path = json_metadata.external_store.root_dir
        metadata_json = base_path / 'metadata.json'
        assert metadata_json.exists()
        with open(metadata_json, mode='r') as f:
            metadata_dict = json.load(f)

        metadata = {key: Metadata(**val)
                    for key, val in metadata_dict.items()}

        assert metadata == json_metadata._metadata_cache
        assert json_metadata['path/to/other.txt'] == meta
        assert len(metadata) == 2

    def test_load_all_metadata(self, json_metadata):
        expected = {'path/to/foo.txt': Metadata(md5='bar')}
        json_metadata._metadata_cache = {}
        loaded = json_metadata.load_all_metadata()
        assert loaded == expected

    def test_load_all_metadata_nofile(self, tmpdir):
        json_metadata = JSONMetadataStore(FileStorage(str(tmpdir)))
        # implicitly called on init anyway
        assert json_metadata._metadata_cache == {}
        # but we also call explicitly
        assert json_metadata.load_all_metadata() == {}

    def test_delete(self, json_metadata):
        assert 'path/to/foo.txt' in json_metadata
        assert len(json_metadata) == 1
        del json_metadata['path/to/foo.txt']
        assert 'path/to/foo.txt' not in json_metadata
        assert len(json_metadata) == 0

    def test_iter(self, json_metadata):
        assert list(json_metadata) == ["path/to/foo.txt"]

    def test_len(self, json_metadata):
        assert len(json_metadata) == 1

    def test_getitem(self, json_metadata):
        assert json_metadata["path/to/foo.txt"] == Metadata(md5="bar")
