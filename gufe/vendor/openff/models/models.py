from collections.abc import Callable
from typing import Any

from openff.units import Quantity
from pydantic import BaseModel

from .types import custom_quantity_encoder


class DefaultModel(BaseModel):
    """A custom Pydantic model used by other components."""

    class Config:
        """Custom Pydantic configuration."""

        json_encoders: dict[Any, Callable] = {
            Quantity: custom_quantity_encoder,
        }
        validate_assignment: bool = True
        arbitrary_types_allowed: bool = True
