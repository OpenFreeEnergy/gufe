# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
import uuid
from typing import List, Optional, Any

from pydantic import BaseModel, PrivateAttr


class ProtocolUnitResult(BaseModel):
    """Result for a single `ProtocolUnit` execution.

    """
    class Config:
        extra = "allow"
        
    _dependencies: List[str] = PrivateAttr()
    _uuid: str = PrivateAttr()

    name: str  # name of the `ProtocolUnit` that produced this `ProtocolUnitResult`
    data: Any  # should likely be fleshed out, currently a free-for-all

    def __init__(self, 
            dependencies: Optional[List["ProtocolUnitResult"]] = None,
            **data):
        """

        """
        dep_uuids = [dep._uuid for dep in dependencies]

        data.update({'_dependencies': dep_uuids})

        super().__init__(**data)
        self._uuid = uuid.uuid4()


class ProtocolDAGResult(BaseModel):
    """Result for a single `ProtocolDAG` execution.

    There may be many of these in a given `ResultStore` for a given `Transformation`.
    Data elements from these objects are combined by `Protocol.gather` into a
    `ProtocolResult`.

    """
    units: List[ProtocolUnitResult]
    
