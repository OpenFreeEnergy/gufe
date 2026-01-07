from pathlib import Path

import pytest

import gufe
from gufe.archival import AlchemicalArchive


def pdr_from_transformation(transformation):
    return gufe.ProtocolDAGResult(protocol_units=[], protocol_unit_results=[], transformation_key=transformation.key)


@pytest.fixture()
def alchemical_archive(benzene_variants_star_map):
    alchemical_network = benzene_variants_star_map
    transformations = sorted(list(alchemical_network.edges))
    # create fake results for the transformations
    transformation_results = {}
    for transformation in transformations:
        transformation_results[transformation.key] = [pdr_from_transformation(transformation)]
    metadata = {"test_meta_key": "test_meta_value"}
    return AlchemicalArchive(
        network=alchemical_network, transformation_results_map=transformation_results, metadata=metadata
    )


class TestArchival:
    def test_roundtrip(self, alchemical_archive):
        json = alchemical_archive.to_json()
        reconstructed = AlchemicalArchive.from_json(content=json)

        assert alchemical_archive == reconstructed

    def test_to_json(self, alchemical_archive, tmpdir):
        # returning a string
        result = alchemical_archive.to_json()
        assert result is not None

        # write to file
        output_file = tmpdir / "archive.json"
        result = alchemical_archive.to_json(output_file)
        assert result is None
        assert output_file.exists()

    def test_from_json(self, alchemical_archive, tmpdir):
        # encode to both a file and a string in memory
        filename = tmpdir / "archive.json"
        alchemical_archive.to_json(filename)
        content = alchemical_archive.to_json()

        assert AlchemicalArchive.from_json(file=filename) == alchemical_archive
        assert AlchemicalArchive.from_json(content=content) == alchemical_archive

        with pytest.raises(ValueError, match="Cannot specify both"):
            AlchemicalArchive.from_json(content=content, file=filename)

    def test_version(self, alchemical_archive):
        # fixture will generate correct version
        assert alchemical_archive.version_gufe == gufe.__version__

    def test_version_warning(self, alchemical_archive):
        alchemical_archive.version_gufe = "0.0.0"
        json = alchemical_archive.to_json()

        with pytest.warns(UserWarning):
            reconstructed = AlchemicalArchive.from_json(content=json)
