# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

"""Long-term storage of ``AlchemicalNetwork`` objects and their results"""

import bisect
import warnings
from typing import Any

import gufe
from gufe.network import AlchemicalNetwork
from gufe.protocols import ProtocolDAGResult
from gufe.tokenization import GufeKey, GufeTokenizable
from gufe.transformations.transformation import TransformationBase


class AlchemicalArchive(GufeTokenizable):
    """A data structure for long-term storage of :class:`.AlchemicalNetwork` objects and the :class:`.ProtocolDAGResult` objects associated with a network's :class:`.Transformation`.

    Parameters
    ----------
    network
        The :class:`.AlchemicalNetwork` being archived.

    transformation_results
        A list of result entries for :class:`.Transformation` objects
        in the provided :class:`.AlchemicalNetwork`. Each entry should
        have the form ``tuple[Transformation,
        list[ProtocolDAGResult]]``.

    metadata
        A dictionary that contains arbitrary metadata associated with
        the archive. This defaults to an empty dictionary.

    version_gufe
        The version of ``gufe`` that produced the archive. By default
        this will be set the version currently installed.
    """

    def __init__(
        self,
        network: AlchemicalNetwork,
        transformation_results: list[tuple[TransformationBase, list[ProtocolDAGResult]]],
        metadata: dict[str, Any] | None = None,
        version_gufe: str | None = None,
    ):
        self._network = network
        self._transformation_results: list[tuple[TransformationBase, list[ProtocolDAGResult]]] = []
        self._validate_transformation_results(transformation_results)

        self._metadata = metadata or {}
        self._version_gufe = version_gufe or gufe.__version__

    @property
    def network(self):
        return self._network

    @property
    def transformation_results(self):
        return self._transformation_results

    @property
    def metadata(self):
        return self._metadata

    @property
    def version_gufe(self):
        return self._version_gufe

    def _validate_transformation_results(self, transformation_results):
        _processed_transformations: set[GufeKey] = set()
        for transformation, pdrs in transformation_results:
            # check that all transformation keys provided are also in
            # the provided AlchemicalNetwork
            if transformation not in self.network.edges:
                raise ValueError(f"{transformation} was not found in {self.network}")

            # only process results if the transformation has not been seen before
            if transformation.key in _processed_transformations:
                msg = f"Duplicate entry for {transformation.key} found in transformation_results"
                raise ValueError(msg)

            _processed_transformations.add(transformation.key)
            entry = [transformation, sorted(set(pdrs))]
            bisect.insort(self._transformation_results, entry, key=lambda e: e[0].key)

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
