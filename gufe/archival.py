import json
import warnings
from dataclasses import asdict, dataclass, field
from os import PathLike
from typing import Any, TextIO

import gufe
from gufe.network import AlchemicalNetwork
from gufe.protocols import ProtocolDAGResult
from gufe.tokenization import JSON_HANDLER, GufeKey, GufeTokenizable, KeyedChain


@dataclass
class AlchemicalArchive:
    network: AlchemicalNetwork
    transformation_results_map: dict[GufeKey, list[ProtocolDAGResult]]
    metadata: dict[str, Any] = field(default_factory=dict)
    version_gufe: str = field(default=gufe.__version__)

    @classmethod
    def from_json(cls, file: PathLike | TextIO | None = None, content: str | None = None):
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
        deserialized["transformation_results_map"] = {
            tokenizable_map[k].key: ProtocolDAGResult.from_keyed_chain(v)
            for k, v in deserialized["transformation_results_map"].items()
        }

        if gufe.__version__ != deserialized["version_gufe"]:
            version_gufe = deserialized["version_gufe"]
            msg = f"Archive was produced with version {version_gufe} and is being loaded with version {gufe.__version__}. If results are unusual, try loading with version {version_gufe}."
            warnings.warn(msg)

        return AlchemicalArchive(**deserialized)

    def to_json(self, file: PathLike | TextIO | None = None) -> None | str:
        dct = asdict(self)

        dct["network"] = dct["network"].to_keyed_chain()
        dct["transformation_results_map"] = {
            k: v.to_keyed_chain() for k, v in dct["transformation_results_map"].items()
        }

        if file is None:
            return json.dumps(dct, cls=JSON_HANDLER.encoder)

        from gufe.utils import ensure_filelike

        with ensure_filelike(file, mode="w") as out:
            json.dump(dct, out, cls=JSON_HANDLER.encoder)

        return None
