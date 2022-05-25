class ResultStore:
    """Class to manage communication between metadata and data storage.

    At this level, we provide an abstraction where client code no longer
    needs to be aware of the nature of the metadata, or even that it exists.
    """
    def __init__(self, external_store, metadata_store):
        self.external_store = external_store
        self.metadata_store = metadata_store

    def store(self, location, byte_data):
        loc, sha2 = self.external_store.create(location, byte_data)
        self.metadata_store.store_metadata(loc, sha2)

    def __iter__(self):
        return iter(self.metadata_store)

    def find_missing_files(self):
        """Identify files listed in metadata but unavailable in storage"""
        return [f for f in self if not self.external_store.exists(f)]

    def load_stream(self, location):
        sha2 = self.metadata_store[location]
        return self.external_store.as_filelike(location, sha2)
