import json
from collections import abc

from typing import Tuple

from gufe.storage.errors import MissingExternalResourceError


class JSONMetadataStore(abc.Mapping):
    # Using JSON for now because it is easy to write this class and doesn't
    # require any external dependencies. It is NOT the right way to go in
    # the long term. API will probably stay the same, though.
    def __init__(self, external_store):
        self.external_store = external_store
        self._metadata_cache = self.load_all_metadata()

    def store_metadata(self, location: str, metadata: str):
        self._metadata_cache[location] = metadata
        metadata_bytes = json.dumps(self._metadata_cache).encode('utf-8')
        _ = self.external_store.store('metadata.json', metadata_bytes)

    def load_all_metadata(self):
        if not self.external_store.exists('metadata.json'):
            return {}

        metadata = self.external_store.get_metadata('metadata.json')
        with self.external_store.load_stream('metadata.json',
                                             metadata) as json_f:
            all_metadata = json.loads(json_f.read().decode('utf-8'))
        return all_metadata

    def __iter__(self):
        return iter(self._metadata_cache)

    def __len__(self):
        return len(self._metadata_cache)

    def __getitem__(self, key):
        return self._metadata_cache[key]
