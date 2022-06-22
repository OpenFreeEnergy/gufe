import json
import abc
import collections

from typing import Tuple

from gufe.storage.errors import MissingExternalResourceError


class MetadataStore(collections.abc.Mapping):
    def __init__(self, external_store):
        self.external_store = external_store
        self._metadata_cache = self.load_all_metadata()

    @abc.abstractmethod
    def store_metadata(self, location: str, metadata: str):
        raise NotImplementedError()

    @abc.abstractmethod
    def load_all_metadata(self):
        raise NotImplementedError()

    @abc.abstractmethod
    def __delitem__(self, location):
        raise NotImplementedError()

    def __getitem__(self, location):
        return self._metadata_cache[location]

    def __iter__(self):
        return iter(self._metadata_cache)

    def __len__(self):
        return len(self._metadata_cache)


class JSONMetadataStore(MetadataStore):
    # Using JSON for now because it is easy to write this class and doesn't
    # require any external dependencies. It is NOT the right way to go in
    # the long term. API will probably stay the same, though.
    def _dump_file(self):
        metadata_bytes = json.dumps(self._metadata_cache).encode('utf-8')
        self.external_store.store_bytes('metadata.json', metadata_bytes)

    def store_metadata(self, location: str, metadata: str):
        self._metadata_cache[location] = metadata
        self._dump_file()

    def load_all_metadata(self):
        if not self.external_store.exists('metadata.json'):
            return {}

        with self.external_store.load_stream('metadata.json') as json_f:
            all_metadata = json.loads(json_f.read().decode('utf-8'))
        return all_metadata

    def __delitem__(self, location):
        del self._metadata_cache[location]
        self._dump_file()
