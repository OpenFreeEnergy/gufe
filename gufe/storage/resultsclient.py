from gufe.storage.resultstore import ResultContainer
from gufe.storage.metadatastore import JSONMetadataStore


class ResultsClient(ResultContainer):
    def __init__(self, external_store):
        # default client is using JSONMetadataStore with the given external
        # result store; users could easily write a subblass that behaves
        # differently
        metadata_store = JSONMetadataStore(external_store)
        self._result_store = _ResultStore(external_store, metadata_store)
        super().__init__(parent=self, path_component=None)

    def store_protocol_dag_result(self, result):
        # I don't know how we get the path information for the protocol dag
        # results
        self.result_store.store(...)

    def _load(self, transformation):
        return TransformationResults(self, transformation)

    # override these two inherited properies since this is always the end of
    # the recursive chain
    @property
    def path(self):
        return 'transformations'

    @property
    def result_store(self):
        return self._result_store


class TransformationResults(_ResultContainer):
    def __init__(self, parent, transformation):
        super().__init__(parent, transformation)
        self.transformation = transformation

    def _load(self, clone):
        return CloneResults(self, clone)


class CloneResults(_ResultContainer):
    def __init__(self, parent, clone):
        super().__init__(parent, clone)
        self.clone = clone

    @staticmethod
    def _to_path_component(item):
        return str(item)

    def _load(self, extension):
        return ExtensionResults(self, extension)


class ExtensionResults(_ResultContainer):
    def __init__(self, parent, extension):
        super().__init__(parent, str(extension))
        self.extension = extension

    @staticmethod
    def _to_path_component(item):
        return str(item)

    def __getitem__(self, filename):
        # different here -- we don't cache the actual file objects
        return self._load(filename)

    def _load(self, filename):
        return self.result_store.load_stream(self.path + "/"
                                             + self.filename)
