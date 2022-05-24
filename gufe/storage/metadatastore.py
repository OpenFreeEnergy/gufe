import json
from collections import abc


class JSONMetadataStore(abc.Mapping):
    # Using JSON for now because it is easy to write this class and doesn't
    # require any external dependencies. It is NOT the right way to go in
    # the long term. API will probably stay the same, though.
    def __init__(self, external_store):
        self.external_store = external_store
        self._metadata_cache = self.load_all_metadata()

    def store_metadata(self, location, sha2):
        self._metadata_cache[location] = sha2
        metadata_bytes = json.dumps(self._metadata_cache).encode('utf-8')
        _ = self.external_store.create('metadata.json', metadata_bytes)

    def load_all_metadata(self):
        if not self.external_store.exists('metadata.json'):
            return {}

        sha2 = self.external_store.get_sha2('metadata.json')
        json_f = self.external_store.as_filelike('metadata.json', sha2)
        metadata = json.loads(json_f.read().decode('utf-8'))
        return metadata

    def __iter__(self):
        return iter(self._metadata_cache)

    def __len__(self):
        return len(self._metadata_cache)

    def __getitem__(self, key):
        try:
            return self._metadata_cache[key]
        except KeyError:
            raise MissingExternalResourceError(f"Hash for '{key}' not "
                                               "found")
