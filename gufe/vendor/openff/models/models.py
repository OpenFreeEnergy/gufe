from collections.abc import Callable
from typing import Any

from openff.units import Quantity
from pydantic.v1 import BaseModel

from .types import custom_quantity_encoder, json_loader


class DefaultModel(BaseModel):
    """A custom Pydantic model used by other components."""

    class Config:
        """Custom Pydantic configuration."""

        json_encoders: dict[Any, Callable] = {
            Quantity: custom_quantity_encoder,
        }
        json_loads: Callable = json_loader
        validate_assignment: bool = True
        arbitrary_types_allowed: bool = True
