import pytest
from unittest import mock

from gufe.storage.externalresource import MemoryStorage
from gufe.storage.resultsclient import (
    ResultsClient, TransformationResults, CloneResults, ExtensionResults
)


@pytest.fixture
def results_client(tmpdir):
    external = MemoryStorage()
    results_client = ResultsClient(external)

    # store one file with contents "foo"
    results_client.result_store.store_bytes(
        "transformations/MAIN_TRANS/0/0/file.txt",
        "foo".encode('utf-8')
    )

    # create some empty files as well
    empty_files = [
        "transformations/MAIN_TRANS/0/0/other.txt",
        "transformations/MAIN_TRANS/0/1/file.txt",
        "transformations/MAIN_TRANS/1/0/file.txt",
        "transformations/OTHER_TRANS/0/0/file.txt",
        "other_dir/file.txt",
    ]

    for file in empty_files:
        results_client.result_store.store_bytes(file, b"")  # empty

    return results_client


def _make_mock_transformation(hash_str):
    return mock.Mock(
        # TODO: fill this in so that it mocks out the digest we use
    )


class _ResultContainerTest:
    @staticmethod
    def get_container(results_client):
        raise NotImplementedError()

    def _getitem_object(self, container):
        raise NotImplementedError()

    def test_iter(self, results_client):
        container = self.get_container(results_client)
        assert set(container) == set(self.expected_files)

    def _get_key(self, as_object, container):
        # TODO: this isn't working yet -- need an interface that allows me
        # to patch the hex digest that we'll be using
        if as_object:
            pytest.skip("Waiting on hex digest patching")
        obj = self._getitem_object(container)
        # next line uses some internal implementation
        key = obj if as_object else obj._path_component
        return key, obj

    @pytest.mark.parametrize('as_object', [True, False])
    def test_getitem(self, as_object, results_client):
        container = self.get_container(results_client)
        key, obj = self._get_key(as_object, container)
        assert container[key] == obj

    @pytest.mark.parametrize('as_object', [True, False])
    def test_div(self, as_object, results_client):
        container = self.get_container(results_client)
        key, obj = self._get_key(as_object, container)
        assert container / key == obj

    @pytest.mark.parametrize('load_with', ['div', 'getitem'])
    def test_caching(self, results_client, load_with):
        # used to test caching regardless of how first loaded was loaded
        container = self.get_container(results_client)
        key, obj = self._get_key(False, container)

        if load_with == 'div':
            loaded = container / key
        elif load_with == 'getitem':
            loaded = container[key]
        else:  # -no-cov-
            raise RuntimeError(f"Bad input: can't load with '{load_with}'")

        assert loaded == obj
        assert loaded is not obj
        reloaded_div = container / key
        reloaded_getitem = container[key]

        assert loaded is reloaded_div
        assert reloaded_div is reloaded_getitem

    def test_load_stream(self, results_client):
        container = self.get_container(results_client)
        loc = "transformations/MAIN_TRANS/0/0/file.txt"
        with container.load_stream(loc) as f:
            assert f.read().decode('utf-8') == "foo"

    def test_load_bytes(self, results_client):
        container = self.get_container(results_client)
        loc = "transformations/MAIN_TRANS/0/0/file.txt"
        assert container.load_bytes(loc).decode('utf-8') == "foo"


    def test_path(self, results_client):
        container = self.get_container(results_client)
        assert container.path == self.expected_path

    def test_result_store(self, results_client):
        container = self.get_container(results_client)
        assert container.result_store == results_client.result_store


class TestResultsClient(_ResultContainerTest):
    expected_files = [
        "transformations/MAIN_TRANS/0/0/file.txt",
        "transformations/MAIN_TRANS/0/0/other.txt",
        "transformations/MAIN_TRANS/0/1/file.txt",
        "transformations/MAIN_TRANS/1/0/file.txt",
        "transformations/OTHER_TRANS/0/0/file.txt",
    ]
    expected_path = "transformations"

    @staticmethod
    def get_container(results_client):
        return results_client

    def _getitem_object(self, container):
        return TransformationResults(
            parent=container,
            transformation=_make_mock_transformation("MAIN_TRANS")
        )

    def test_store_protocol_dag_result(self):
        pytest.skip("Not implemented yet")

    def test_delete(self, results_client):
        file_to_delete = self.expected_files[0]
        storage = results_client.result_store.external_store
        assert storage.exists(file_to_delete)
        results_client.delete(file_to_delete)
        assert not storage.exists(file_to_delete)


class TestTransformationResults(_ResultContainerTest):
    expected_files = [
        "transformations/MAIN_TRANS/0/0/file.txt",
        "transformations/MAIN_TRANS/0/0/other.txt",
        "transformations/MAIN_TRANS/0/1/file.txt",
        "transformations/MAIN_TRANS/1/0/file.txt",
    ]
    expected_path = "transformations/MAIN_TRANS"

    @staticmethod
    def get_container(results_client):
        container = TransformationResults(
            parent=TestResultsClient.get_container(results_client),
            transformation=_make_mock_transformation("MAIN_TRANS")
        )
        container._path_component = "MAIN_TRANS"
        return container

    def _getitem_object(self, container):
        return CloneResults(parent=container, clone=0)


class TestCloneResults(_ResultContainerTest):
    expected_files = [
        "transformations/MAIN_TRANS/0/0/file.txt",
        "transformations/MAIN_TRANS/0/0/other.txt",
        "transformations/MAIN_TRANS/0/1/file.txt",
    ]
    expected_path = "transformations/MAIN_TRANS/0"

    @staticmethod
    def get_container(results_client):
        return CloneResults(
            parent=TestTransformationResults.get_container(results_client),
            clone=0
        )

    def _getitem_object(self, container):
        return ExtensionResults(parent=container, extension=0)


class TestExtensionResults(_ResultContainerTest):
    expected_files = [
        "transformations/MAIN_TRANS/0/0/file.txt",
        "transformations/MAIN_TRANS/0/0/other.txt",
    ]
    expected_path = "transformations/MAIN_TRANS/0/0"

    @staticmethod
    def get_container(results_client):
        return ExtensionResults(
            parent=TestCloneResults.get_container(results_client),
            extension=0
        )

    def _get_key(self, as_object, container):
        if self.as_object:  # -no-cov-
            raise RuntimeError("TestExtensionResults does not support "
                               "as_object=True")
        path = "transformations/MAIN_TRANS/0/0/"
        fname = "file.txt"
        return fname, container.result_store.load_stream(path + fname)

    # things involving div and getitem need custom treatment
    def test_div(self, results_client):
        container = self.get_container(results_client)
        with container / "file.txt" as f:
            assert f.read().decode('utf-8') == "foo"

    def test_getitem(self, results_client):
        container = self.get_container(results_client)
        with container["file.txt"] as f:
            assert f.read().decode('utf-8') == "foo"

    def test_caching(self, results_client):
        # this one does not cache results; the cache should remain empty
        container = self.get_container(results_client)
        assert container._cache == {}
        from_div = container / "file.txt"
        assert container._cache == {}
        from_getitem = container["file.txt"]
        assert container._cache == {}
