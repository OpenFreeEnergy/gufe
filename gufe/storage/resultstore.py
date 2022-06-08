from gufe.storage.errors import MissingExternalResourceError


class ResultStore:
    """Class to manage communication between metadata and data storage.

    At this level, we provide an abstraction where client code no longer
    needs to be aware of the nature of the metadata, or even that it exists.
    """
    def __init__(self, external_store, metadata_store):
        self.external_store = external_store
        self.metadata_store = metadata_store

    def store(self, location, byte_data):
        _, metadata = self.external_store.store(location, byte_data)
        self.metadata_store.store_metadata(location, metadata)

    def delete(self, location):
        del self.metadata_store[location]
        self.external_store.delete(location)

    def __iter__(self):
        return iter(self.metadata_store)

    def find_missing_files(self):
        """Identify files listed in metadata but unavailable in storage"""
        return [f for f in self if not self.external_store.exists(f)]

    def load_stream(self, location):
        try:
            metadata = self.metadata_store[location]
        except KeyError:
            raise MissingExternalResourceError(f"Hash for '{location}' not "
                                               "found")

        return self.external_store.load_stream(location, metadata)
