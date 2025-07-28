# Vendored from https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/pydantic.py
"""Pydantic base model with custom settings."""

from typing import Any

from pydantic import BaseModel, ConfigDict


class _BaseModel(BaseModel):
    """A custom Pydantic model used by other components."""

    model_config = ConfigDict(
        validate_assignment=True,
        arbitrary_types_allowed=True,
    )

    def model_dump(self, **kwargs) -> dict[str, Any]:
        return super().model_dump(serialize_as_any=True, **kwargs)

    def model_dump_json(self, **kwargs) -> str:
        return super().model_dump_json(serialize_as_any=True, **kwargs)

    def __eq__(self, other: Any) -> bool:
        # reproduces pydantic v1 equality, since v2 checks for private attr equality,
        # which results in frozen/unfrozen objects not being equal
        # https://github.com/pydantic/pydantic/blob/2486e068e85c51728c9f2d344cfee2f7e11d555c/pydantic/v1/main.py#L911
        if isinstance(other, BaseModel):
            return self.model_dump() == other.model_dump()
        else:
            return self.model_dump() == other
