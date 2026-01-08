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
            dct[key] = tuple(sorted(value))


@dataclass
class AlchemicalArchive:
    network: AlchemicalNetwork
    transformation_results_map: dict[GufeKey, tuple[ProtocolDAGResult]]
    metadata: dict[str, Any] = field(default_factory=dict)
    version_gufe: str = field(default=gufe.__version__)

    def __post_init__(self):
        # ensure that ProtocolDAGResults are always in order by gufe key
        _dict_sort_values(self.transformation_results_map)

    @classmethod
    def from_json(cls, file: PathLike | TextIO | None = None, content: str | None = None):
        """Create an ``AlchemicalArchive`` from JSON.

        Parameters
        ----------
        file
            A filepath or file-like object to read JSON from.

        content
            JSON content as a string.

        Returns
        -------
        AlchemicalArchive
        """
        if content is not None and file is not None:
            raise ValueError("Cannot specify both `content` and `file`; only one input allowed")
        elif content is None and file is None:
            raise ValueError("Must specify either `content` and `file` for JSON input")

        if content is not None:
            deserialized = json.loads(content, cls=JSON_HANDLER.decoder)
        else:
            from gufe.utils import ensure_filelike

            with ensure_filelike(file, mode="r") as f:
                deserialized = json.load(f, cls=JSON_HANDLER.decoder)

        # provide our own tokenization map for later updating possibly
        # stale transformation keys
        tokenizable_map: dict[str, GufeTokenizable] = {}
        deserialized["network"] = KeyedChain(deserialized["network"]).to_gufe(tokenizable_map=tokenizable_map)
        for key, results in deserialized["transformation_results_map"].items():
            deserialized["transformation_results_map"][tokenizable_map[key].key] = [
                ProtocolDAGResult.from_keyed_chain(v) for v in results
            ]

        if gufe.__version__ != deserialized["version_gufe"]:
            version_gufe = deserialized["version_gufe"]
            msg = f"Archive was produced with version {version_gufe} and is being loaded with version {gufe.__version__}. If results are unusual, try loading with version {version_gufe}."
            warnings.warn(msg)

        return AlchemicalArchive(**deserialized)

    def to_json(self, file: PathLike | TextIO | None = None) -> None | str:
        """Encode the ``AlchemicalArchive`` to JSON.

        Parameters
        ----------
        file
            A filepath or file-like object to write JSON to.

        Returns
        -------
        str | None
            The JSON encoded ``AlchemicalArchive`` if ``file`` was provided; otherwise ``None``.
        """

        dct = asdict(self)

        dct["network"] = dct["network"].to_keyed_chain()
        for key, results in dct["transformation_results_map"].items():
            dct["transformation_results_map"][key] = [result.to_keyed_chain() for result in results]

        if file is None:
            return json.dumps(dct, sort_keys=True, cls=JSON_HANDLER.encoder)

        from gufe.utils import ensure_filelike

        with ensure_filelike(file, mode="w") as out:
            json.dump(dct, out, sort_keys=True, cls=JSON_HANDLER.encoder)

        return None

    @cached_property
    def md5(self):
        dumped = self.to_json()
        hasher = hashlib.md5(dumped.encode(), usedforsecurity=False)
        return hasher.hexdigest()

    def __hash__(self):
        return hash(self.md5)
