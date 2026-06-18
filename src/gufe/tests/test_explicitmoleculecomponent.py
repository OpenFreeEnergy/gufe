import pickle
from unittest import mock

import pytest


def _pickle_roundtrip(obj):
    return pickle.loads(pickle.dumps(obj))


def _dict_roundtrip(obj):
    return type(obj).from_dict(obj.to_dict())


def _msgpack_roundtrip(obj):
    return type(obj).from_msgpack(content=obj.to_msgpack())


class ExplicitMoleculeComponentMixin:
    @pytest.mark.parametrize(
        "roundtrip",
        [
            pytest.param(_pickle_roundtrip, id="pickle"),
            pytest.param(_dict_roundtrip, id="dict"),
            pytest.param(_msgpack_roundtrip, id="msgpack"),
        ],
    )
    def test_name_and_equality_after_round_trip(self, instance, roundtrip):
        new_instance = roundtrip(instance)

        assert new_instance == instance
        assert new_instance.key == instance.key
        assert type(new_instance) is type(instance)
        assert new_instance is instance

        assert new_instance.name == instance.name

        mol = new_instance.to_rdkit()
        assert mol.HasProp("ofe-name")
        assert mol.GetProp("ofe-name") == instance.name

    def test_pickle_roundtrip_reconstructs_molecule_props_with_empty_registry(self, instance):
        payload = pickle.dumps(instance)

        patch_loc = "gufe.tokenization.TOKENIZABLE_REGISTRY"
        with mock.patch.dict(patch_loc, {}, clear=True):
            new_instance = pickle.loads(payload)

        assert new_instance == instance
        assert new_instance is not instance

        mol = new_instance.to_rdkit()
        assert mol.HasProp("ofe-name")
        assert mol.GetProp("ofe-name") == instance.name
