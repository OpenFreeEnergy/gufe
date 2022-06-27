# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import abc
from typing import Any

from .resultserver import ResultServer
from .metadatastore import JSONMetadataStore


class _ResultContainer(abc.ABC):
    """
    Abstract class, represents all data under some level of the heirarchy.
    """
    def __init__(self, parent, path_component):
        self.parent = parent
        self._path_component = self._to_path_component(path_component)
        self._cache = {}

    def __eq__(self, other):
        return (
            isinstance(other, self.__class__)
            and self.path == other.path
        )

    @staticmethod
    def _to_path_component(item: Any) -> str:
        """Convert input (object or string) to path string"""
        if isinstance(item, str):
            return item

        # TODO: instead of str(hash(...)), this should return the digest
        # that is being introduced in another PR; Python hash is not stable
        # across sessions
        return str(hash(item))

    def __getitem__(self, item):
        # code for the case this is a file
        if item in self.result_server:
            return self.result_server.load_stream(item)

        # code for the case this is a "directory"
        hash_item = self._to_path_component(item)

        if hash_item not in self._cache:
            self._cache[hash_item] = self._load_next_level(item)

        return self._cache[hash_item]

    def __truediv__(self, item):
        return self[item]

    @abc.abstractmethod
    def _load_next_level(self, item):
        raise NotImplementedError()

    def __iter__(self):
        for loc in self.result_server:
            if loc.startswith(self.path):
                yield loc

    def load_stream(self, location, *, allow_changed=False):
        return self.result_server.load_stream(location, allow_changed)

    def load_bytes(self, location, *, allow_changed=False):
        with self.load_stream(location, allow_changed=allow_changed) as f:
            byte_data = f.read()

        return byte_data

    @property
    def path(self):
        return self.parent.path + "/" + self._path_component

    @property
    def result_server(self):
        return self.parent.result_server

    def __repr__(self):
        # probably should include repr of external store, too
        return f"{self.__class__.__name__}({self.path})"


class ResultClient(_ResultContainer):
    def __init__(self, external_store):
        # default client is using JSONMetadataStore with the given external
        # result store; users could easily write a subblass that behaves
        # differently
        metadata_store = JSONMetadataStore(external_store)
        self._result_server = ResultServer(external_store, metadata_store)
        super().__init__(parent=self, path_component=None)

    def delete(self, location):
        self._result_server.delete(location)

    def store_protocol_dag_result(self, result):
        # I don't know how we get the path information for the protocol dag
        # results
        self.result_server.store(...)

    def _load_next_level(self, transformation):
        return TransformationResult(self, transformation)

    # override these two inherited properies since this is always the end of
    # the recursive chain
    @property
    def path(self):
        return 'transformations'

    @property
    def result_server(self):
        return self._result_server


class TransformationResult(_ResultContainer):
    def __init__(self, parent, transformation):
        super().__init__(parent, transformation)
        self.transformation = transformation

    def _load_next_level(self, clone):
        return CloneResult(self, clone)


class CloneResult(_ResultContainer):
    def __init__(self, parent, clone):
        super().__init__(parent, clone)
        self.clone = clone

    @staticmethod
    def _to_path_component(item):
        return str(item)

    def _load_next_level(self, extension):
        return ExtensionResult(self, extension)


class ExtensionResult(_ResultContainer):
    def __init__(self, parent, extension):
        super().__init__(parent, str(extension))
        self.extension = extension

    @staticmethod
    def _to_path_component(item):
        return str(item)

    def __getitem__(self, filename):
        # different here -- we don't cache the actual file objects
        return self._load_next_level(filename)

    def _load_next_level(self, filename):
        return self.result_server.load_stream(self.path + "/" + filename)
