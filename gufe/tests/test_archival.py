import hashlib
from pathlib import Path

import pytest

import gufe
from gufe.archival import AlchemicalArchive

from .test_tokenization import GufeTokenizableTestsMixin


def pdr_from_transformation(transformation):
    return gufe.ProtocolDAGResult(protocol_units=[], protocol_unit_results=[], transformation_key=transformation.key)


class TestArchival(GufeTokenizableTestsMixin):

    @pytest.fixture()
    def instance(benzene_variants_star_map):
        alchemical_network = benzene_variants_star_map
        transformations = sorted(list(alchemical_network.edges))
        # create fake results for the transformations
        transformation_results = {}
        for transformation in transformations:
            transformation_results[transformation.key] = [pdr_from_transformation(transformation)]
        metadata = {"test_meta_key": "test_meta_value", "meta_ordered": [3, 2, 1]}
        return AlchemicalArchive(
            network=alchemical_network, transformation_results_map=transformation_results, metadata=metadata
        )

    def test_roundtrip(self, instance):
        json = instance.to_json()
        reconstructed = AlchemicalArchive.from_json(content=json)

        assert instance == reconstructed

    def test_to_json(self, instance, tmpdir):
        # returning a string
        result = instance.to_json()
        assert result is not None

        # write to file
        output_file = tmpdir / "archive.json"
        result = instance.to_json(output_file)
        assert result is None
        assert output_file.exists()

    def test_from_json(self, instance, tmpdir):
        # encode to both a file and a string in memory
        filename = tmpdir / "archive.json"
        instance.to_json(filename)
        content = instance.to_json()

        assert AlchemicalArchive.from_json(file=filename) == instance
        assert AlchemicalArchive.from_json(content=content) == instance

        with pytest.raises(ValueError, match="Cannot specify both"):
            AlchemicalArchive.from_json(content=content, file=filename)

    def test_version(self, instance):
        # fixture will generate correct version
        assert instance.version_gufe == gufe.__version__

    def test_version_warning(self, instance):
        instance.version_gufe = "0.0.0"
        json = instance.to_json()

        with pytest.warns(UserWarning):
            reconstructed = AlchemicalArchive.from_json(content=json)

    def test_metadata(self, instance):
        assert instance.metadata == {"test_meta_key": "test_meta_value", "meta_ordered": [3, 2, 1]}
