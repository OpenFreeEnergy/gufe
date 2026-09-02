import importlib
import warnings
from unittest import mock

import pytest

import gufe
from gufe.archival import AlchemicalArchive

from .test_tokenization import GufeTokenizableTestsMixin


def pdr_from_transformation(transformation):
    return gufe.ProtocolDAGResult(protocol_units=[], protocol_unit_results=[], transformation_key=transformation.key)


class TestArchival(GufeTokenizableTestsMixin):
    cls = AlchemicalArchive
    repr = None

    @pytest.fixture()
    def instance(instance, alchem_network_benzene_variants):
        alchemical_network = alchem_network_benzene_variants
        transformations = sorted(list(alchemical_network.edges))
        # create fake results for the transformations
        transformation_results = []
        for transformation in transformations:
            transformation_results.append([transformation, [pdr_from_transformation(transformation)]])
        metadata = {"test_meta_key": "test_meta_value", "meta_ordered": [3, 2, 1]}
        instance = AlchemicalArchive(
            network=alchemical_network, transformation_results=transformation_results, metadata=metadata
        )
        return instance

    def test_version(self, instance):
        # fixture will generate correct version
        assert instance.version_gufe == gufe.__version__

    # only warn for a difference in major or minor version,
    # do not warn for difference in patch or dev version
    @pytest.mark.parametrize(
        "archive_ver,current_ver,should_warn",
        [
            ("1.12.0", "1.13.0", True),
            ("1.12.0", "1.12.1", False),
            ("1.12.0a0", "1.12.0", False),
            ("1.12.0-a0", "1.12.0", False),
            ("1.12.1", "1.13", True),
            ("1.13", "2.0", True),
        ],
    )
    def test_version_warning(
        self,
        instance,
        archive_ver,
        current_ver,
        should_warn,
    ):

        def json_with_archive_ver(new_version):
            new_instance = instance.copy_with_replacements(version_gufe=new_version)
            return new_instance.to_json()

        # pretend this json was created with `archive_ver`
        json = json_with_archive_ver(archive_ver)
        with mock.patch("gufe.__version__", current_ver):
            if should_warn:
                with pytest.warns(UserWarning, match="try loading with"):
                    _ = AlchemicalArchive.from_json(content=json)
            else:
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter("always")
                    _ = AlchemicalArchive.from_json(content=json)
                    assert len(w) == 0

    def test_metadata(self, instance):
        assert instance.metadata == {"test_meta_key": "test_meta_value", "meta_ordered": [3, 2, 1]}

    def test_invalid_transformation_key(self, instance):
        valid_transformations = sorted(list(instance.network.edges))
        invalid_transformation = valid_transformations[0].copy_with_replacements(name="invalid_transformation")
        invalid_pdr = pdr_from_transformation(invalid_transformation)

        with pytest.raises(ValueError, match=r"^.+ was not found in"):
            instance.copy_with_replacements(
                transformation_results=instance.transformation_results + [[invalid_transformation, [invalid_pdr]]]
            )

    def test_repeated_transformation(self, instance):
        new_results = instance.transformation_results.copy()
        new_results += [new_results[-1]]

        with pytest.raises(ValueError, match="Duplicate entry for"):
            instance.copy_with_replacements(transformation_results=new_results)

    def test_transformation_ordering(self, instance):
        new_results = instance.transformation_results[::-1]
        reconstructed = instance.copy_with_replacements(transformation_results=new_results)

        assert new_results != instance.transformation_results
        assert instance == reconstructed

    def test_immutability(self, instance):
        with pytest.raises(AttributeError, match="has no setter"):
            instance.version_gufe = "25"


def test_regression_archive_serialization():
    with importlib.resources.path("gufe.tests.data", "alchemical_archive.json") as file:
        filename = str(file)

    archive = AlchemicalArchive.from_json(file=filename)

    assert archive.metadata == {"test_meta_key": "test_meta_value", "meta_ordered": [3, 2, 1]}
    assert len(archive.network.edges) == len(archive.transformation_results) == 12

    for transformation, pdrs in archive.transformation_results:
        assert transformation in archive.network.edges
        assert len(pdrs) == 1

    assert archive.version_gufe == "1.7.1.dev46+gb75e1476f.d20260203"
