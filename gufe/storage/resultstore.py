class ResultStore:
    """Facade class to abstract out saving/loading data.
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


class _ResultContainer:
    """
    Abstract class, represents all data under some level of the heirarchy.
    """
    def __init__(self, parent, path_component):
        self.parent = parent
        self._path_component = self._to_path_component(path_component)
        self._cache = {}

    @staticmethod
    def _to_path_component(item):
        if isinstance(item, str):
            return item
        return str(hash(item))

    def __getitem__(self, item):
        hash_item = self._to_path_component(item)

        if hash_item not in self._cache:
            self._cache[hash_item] = self._load(item)

        return self._cache[hash_item]

    def __truediv__(self, item):
        return self[item]

    def _load(self, item):
        raise NotImplementedError()

    def __iter__(self):
        for loc in self.result_store:
            if loc.startswith(self.path):
                yield loc

    def load_stream(self, location):
        return self.result_store.load_stream(location)

    @property
    def path(self):
        return self.parent.path + "/" + self._path_component

    @property
    def result_store(self):
        return self.parent.result_store

    def __repr__(self):
        # probably should include repr of external store, too
        return f"{self.__class__.__name__}({self.path})"


