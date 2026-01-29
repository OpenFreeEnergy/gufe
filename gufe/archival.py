import warnings
from typing import Any

import gufe
from gufe.network import AlchemicalNetwork
from gufe.protocols import ProtocolDAGResult
from gufe.tokenization import GufeKey, GufeTokenizable
from gufe.transformations.transformation import TransformationBase


class AlchemicalArchive(GufeTokenizable):
    def __init__(
        self,
        network: AlchemicalNetwork,
        transformation_results: list[tuple[TransformationBase, list[ProtocolDAGResult]]],
        metadata: dict[str, Any],
        version_gufe: str | None = None,
    ):
        self.network = network

        network_transformation_keys = {edge.key for edge in self.network.edges}

        _results: dict[TransformationBase, list[ProtocolDAGResult]] = {}
        for transformation, pdrs in transformation_results:
            # check that all transformation keys provided are also in
            # the provided AlchemicalNetwork
            if transformation.key not in network_transformation_keys:
                raise ValueError(f"{transformation} was not found in {self.network}")

            # TODO: this needs more testing
            if previous := _results.get(transformation):
                _results[transformation] = sorted(set(previous + list(pdrs)))
            else:
                _results[transformation] = sorted(set(pdrs))

        self.transformation_results = []
        for transformation in sorted(_results.keys()):
            self.transformation_results.append([transformation, _results[transformation]])

        self.metadata = metadata
        self.version_gufe = version_gufe or gufe.__version__

    def _to_dict(self):
        return {
            "network": self.network,
            "transformation_results": self.transformation_results,
            "metadata": self.metadata,
            "version_gufe": self.version_gufe,
        }

    @classmethod
    def _defaults(cls):
        return {}

    @classmethod
    def _from_dict(cls, dct):
        returned_instance = cls(**dct)
        if gufe.__version__ != returned_instance.version_gufe:
            version_gufe = returned_instance.version_gufe
            msg = f"Archive was produced with version {version_gufe} and is being loaded with version {gufe.__version__}. If results are unusual, try loading with version {version_gufe}."
            warnings.warn(msg)

        return returned_instance
