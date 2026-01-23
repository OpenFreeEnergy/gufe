import warnings
from typing import Any

import gufe
from gufe.network import AlchemicalNetwork
from gufe.protocols import ProtocolDAGResult
from gufe.tokenization import GufeKey, GufeTokenizable


class AlchemicalArchive(GufeTokenizable):
    def __init__(
        self,
        network: AlchemicalNetwork,
        transformation_results_map: dict[GufeKey | str, list[ProtocolDAGResult]],
        metadata: dict[str, Any],
        version_gufe: str | None = None,
    ):
        self.network = network
        self.transformation_results_map = {}

        network_transformation_keys = [str(edge.key) for edge in self.network.edges]

        for transformation, results in transformation_results_map.items():
            # check that all keys in transformation map are GufeKeys
            if not isinstance(transformation, str):
                raise ValueError(f"Keys of transformation_results_map must be instances of GufeKey or str")

            # check that all transformation keys provided are also in
            # the provided AlchemicalNetwork
            if transformation not in network_transformation_keys:
                raise ValueError(f"{transformation} was not found in {self.network}")

            self.transformation_results_map[str(transformation)] = sorted(list(results))

        self.metadata = metadata
        self.version_gufe = version_gufe or gufe.__version__

    def _to_dict(self):
        return {
            "network": self.network,
            "transformation_results_map": self.transformation_results_map,
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
