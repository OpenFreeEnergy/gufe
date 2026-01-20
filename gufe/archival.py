import hashlib
import json
import warnings
from dataclasses import asdict, dataclass, field
from functools import cached_property
from os import PathLike
from typing import Any, TextIO

import gufe
from gufe.network import AlchemicalNetwork
from gufe.protocols import ProtocolDAGResult
from gufe.tokenization import JSON_HANDLER, GufeKey, GufeTokenizable, KeyedChain


def _dict_sort_values(dct):
    for key, value in dct.items():
        if isinstance(value, (list, tuple)):
            dct[key] = list(sorted(value))


class AlchemicalArchive(GufeTokenizable):
    def __init__(self, network, transformation_results_map, metadata, version_gufe=None):
        self.network = network
        self.transformation_results_map = transformation_results_map
        _dict_sort_values(self.transformation_results_map)
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
