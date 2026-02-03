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
    def instance(instance, benzene_variants_star_map):
        alchemical_network = benzene_variants_star_map
        transformations = sorted(list(alchemical_network.edges))
        # create fake results for the transformations
        transformation_results = []
        for transformation in transformations:
            transformation_results.append([transformation, [pdr_from_transformation(transformation)]])
        metadata = {"test_meta_key": "test_meta_value", "meta_ordered": [3, 2, 1]}
        return AlchemicalArchive(
            network=alchemical_network, transformation_results=transformation_results, metadata=metadata
        )

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
