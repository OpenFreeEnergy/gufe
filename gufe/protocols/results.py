# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from typing import List

from pydantic import BaseModel


class ProtocolUnitResult(BaseModel):
    """Result for a single `ProtocolUnit` execution.

    """
    class Config:
        extra = "allow"


class ProtocolDAGResult(BaseModel):
    """Result for a single `ProtocolDAG` execution.

    There may be many of these in a given `ResultStore` for a given `Transformation.protocol`.

    """
    units: List[ProtocolUnitResult]
    


class ProtocolResult(BaseModel):
    """Container for all `ProtocolDAGResult`s for a given `Transformation.protocol`.

    There will be exactly one of these in a given `ResultStore` for a given `Transformation.protocol`.

    """
    class Config:
        extra = "allow"
