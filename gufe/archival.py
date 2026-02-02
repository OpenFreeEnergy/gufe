import bisect
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
        self.transformation_results: list[tuple[TransformationBase, list[ProtocolDAGResult]]] = []
        self._validate_transformation_results(transformation_results)

        self.metadata = metadata
        self.version_gufe = version_gufe or gufe.__version__

    def _validate_transformation_results(self, transformation_results):
        _processed_transformations: set[GufeKey] = set()
        for transformation, pdrs in transformation_results:
            # check that all transformation keys provided are also in
            # the provided AlchemicalNetwork
            if transformation not in self.network.edges:
                raise ValueError(f"{transformation} was not found in {self.network}")

            # only process results if the transformation has not been seen before
            if transformation in _processed_transformations:
                msg = f"Duplicate entry for {transformation.key} found in transformation_results"
                raise ValueError(msg)

            _processed_transformations.add(transformation.key)
            entry = [transformation, sorted(set(pdrs))]
            bisect.insort(self.transformation_results, entry, key=lambda e: e[0].key)

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
